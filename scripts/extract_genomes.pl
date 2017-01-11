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
use Shrub::GTO;
use File::Copy::Recursive;

=head1 Extract Genome Data to a Directory

    extract_genomes.pl [ options ] outDir

This script accepts a list of genome IDs as input and extracts FASTA or GTO information to an output directory.
The genome files created have the same name as the genome ID with a suffix of C<.fa> or C<.gto> accordingly.
The directories thus created can be used for processing by mass quality-check scripts, among other things.

=head2 Parameters

The single positional parameter is the name of the output directory.

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item fasta

If specified, FASTA files will be output. The default is to output GTO files.

=item missing

If specified, only files not already present in the output directory will be created.

=item col

The index (1-based) of the input column containing the genome IDs. The default is C<0>, indicating the
last column.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['fasta', 'create FASTA files'],
        ['missing', 'process new genomes only'],
        ['col|c=i', 'input column (1-based)', { default => 0 }]
        );
# Create a statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir.";
}
# Get the options.
my $fastaFormat = $opt->fasta;
my $missing = $opt->missing;
my $col = $opt->col - 1;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the genome IDs.
while (! eof $ih) {
    # Read the next genome ID.
    my $line = <$ih>;
    $line =~ s/[\r\n]+//g;
    $stats->Add(lineIn => 1);
    my @fields = split /\t/, $line;
    my $genomeID = $fields[$col];
    # Compute the output file name.
    my $outFileName = "$outDir/$genomeID." . ($fastaFormat ? 'fa' : 'gto');
    # If the file exists and this is missing-mode, skip it.
    if (-s $outFileName && $missing) {
        $stats->Add(skipFound => 1);
        print "$genomeID already in directory.\n";
    } elsif ($fastaFormat) {
        # Here we are producing a FASTA file. We need to ask for the repo file.
        my $fromFile = $shrub->genome_fasta($genomeID);
        if (! $fromFile) {
            print "$genomeID not in database.\n";
            $stats->Add(skipNoSuchGenome => 1);
        } else {
            print "Copying $fromFile to $outFileName.\n";
            File::Copy::Recursive::fcopy($fromFile, $outFileName) || die "Could not copy $fromFile: $!";
            $stats->Add(fastaCopied => 1);
        }
    } else {
        # Here we are producing a GTO file.
        my $gto = Shrub::GTO->new($shrub, $genomeID);
        if (! $gto) {
            print "$genomeID not in database.\n";
            $stats->Add(skipNoSuchGenome => 1);
        } else {
            print "Storing $outFileName.\n";
            $gto->destroy_to_file($outFileName);
            $stats->Add(gtoStored => 1);
        }
    }
}
print "All done.\n" . $stats->Show();