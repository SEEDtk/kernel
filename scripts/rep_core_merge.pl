#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


use strict;
use warnings;
use FIG_Config;
use Shrub;
use ScriptUtils;
use Stats;
use File::Copy::Recursive;

=head1 merge Core Genomes into Representative Genome Files

    rep_core_merge.pl [ options ] outDir

This script will examine the genome list in a representative genome directory and add the quality core genomes
to the two files. The C<complete.genomes> file is read to get a list of the genomes already processed, then
the Shrub database is read to get a list of the core genomes not already in the list. These genomes are
appended to the C<complete.genomes> file and then their seed protein sequences are added to the
C<6.1.1.20.fasta> file.

=head2 Parameters

The positional parameter is the name of the directory containing the representative genome files.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item minlen

The minimum acceptable length for the seed protein. The default is 209.

=item maxlen

The maximum acceptable length for the seed protein. The default is 485.


=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ['minlen=i', 'minimum protein length', { default => 209 }],
        ['maxlen=i', 'maximum protein length', { default => 485 }],
        );
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the genome file directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Invalid or missing output directory $outDir.";
}
# Get the length limits.
my $minlen = $opt->minlen;
my $maxlen = $opt->maxlen;
# Create a hash of the current genomes.
my %current;
my $genomeFile = "$outDir/complete.genomes";
if (! -f $genomeFile) {
    print "No current genome file.\n";
} else {
    open(my $ih, '<', $genomeFile) || die "Could not open $genomeFile: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(genomeLineIn => 1);
        if ($line =~ /^(\d+\.\d+)/) {
            $current{$1} = 1;
            $stats->Add(genomeIn => 1);
        }
    }
    print scalar(keys %current) . " genomes already in file.\n";
    close $ih;
}
# Copy the old files.
print "Copying old files.\n";
my $fastaFile = "$outDir/6.1.1.20.fasta";
File::Copy::Recursive::fcopy($genomeFile, "$genomeFile.new") || die "Could not copy $genomeFile: $!";
File::Copy::Recursive::fcopy($fastaFile, "$fastaFile.new") || die "Could not copy $fastaFile: $!";
# Open the files for appending.
open(my $gh, '>>', "$genomeFile.new") || die "Could not open genome file for output: $!";
open(my $fh, '>>', "$fastaFile.new") || die "Could not open fasta file for output: $!";
# Loop through the good core genomes, collecting the proteins.
print "Scanning the database.\n";
my $q = $shrub->Get('Genome Genome2Feature Feature Feature2Function AND Feature Protein',
        'Genome(well-behaved) = ? AND Feature2Function(security) = ? AND Feature2Function(to-link) = ?',
        [1, 2, 'PhenTrnaSyntAlph'], 'id name Feature(id) Protein(sequence)');
my %gData;
while (my $fData = $q->Fetch()) {
    my ($genome, $name, $fid, $prot) = $fData->Values('id name Feature(id) Protein(sequence)');
    $stats->Add(protInDatabase => 1);
    if (! $current{$genome}) {
        $stats->Add(protInNewGenome => 1);
        push @{$gData{$genome}}, [$fid, $name, $prot];
    }
}
# Now loop through the proteins, writing the good ones.
print "Processing the proteins.\n";
for my $genome (sort keys %gData) {
    # Get the list of protein triples.
    my $prots = $gData{$genome};
    if (@$prots > 1) {
        # Reject if more than one copy of the protein.
        $stats->Add(protTooMany => 1);
    } else {
        # Get the protein's data.
        my $prot = $prots->[0];
        my ($fid, $name, $seq) = @$prot;
        # Verify the length.
        my $len = length $seq;
        if ($len < $minlen) {
            $stats->Add(protTooShort => 1);
        } elsif ($len > $maxlen) {
            $stats->Add(protTooLong => 1);
        } else {
            # Everything is good. Write it out.
            print $fh ">$genome $fid\t$name\n$seq\n";
            print $gh "$genome\t$name\n";
            $stats->Add(protKept => 1);
        }
    }
}
# Copy the files back.
print "Restoring file names.\n";
close $fh; close $gh;
File::Copy::Recursive::fmove("$genomeFile.new", $genomeFile) || die "Could not restore $genomeFile: $!";
File::Copy::Recursive::fmove("$fastaFile.new", $fastaFile) || die "Could not restore $fastaFile: $!";
# All done.
print "All done:\n" . $stats->Show();