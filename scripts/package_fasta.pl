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
use ScriptUtils;
use Stats;
use File::Copy::Recursive;
use RASTlib;
use GenomeTypeObject;
use Shrub;

=head1 Create Genome Packages Out of FASTA Files

    package_fasta.pl [ options ] inDir outDir

This script processes all the genome FASTA files in a directory, and produces genome packages from them.
The input FASTA file should have as its name a genome ID followed by the extension C<.fa>. The genome
package created will be marked as derived from that genome. To create the package, the FASTA file will
be run through RAST and then a C<data.tbl> file created.

=head2 Parameters

The positional parameters are the names of the input directory containing the FASTA files and the output
directory to contain the genome packages created.

The command-line options are those given in L<Shrub::script_options> plus the following.

=over 4

=over 4

=item user

User name for RAST access.

=item password

Password for RAST access.

=item sleep

Sleep interval in seconds while waiting for the job to complete. The default is C<60>.

=back

=cut

$| = 1;
my $stats = Stats->new();
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outDir', Shrub::script_options(),
        ["user|u=s", "user name for RAST access", { default => $ENV{RASTUSER} }],
        ["password|p=s", "password for RAST access", { default => $ENV{RASTPASS} }],
        ["sleep=i", "sleep interval for status polling", { default => 60 }],
        );
# Get the positional parameters.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Invalid or missing output directory $outDir.";
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get all the FASTA files.
opendir(my $dh, $inDir) || die "Could not open input directory: $!";
my @files = grep { $_ =~ /^\d+\.\d+\.fa$/ } readdir $dh;
closedir $dh;
# Loop through them.
for my $file (@files) {
    # Get the source info.
    my ($genomeID, $taxonID) = ($file =~ /^((\d+)\.\d+)/);
    print "Processing $file.\n";
    # Open the FASTA file.
    $stats->Add(fastaFileIn => 1);
    open(my $ih, '<', "$inDir/$file") || die "Could not open $file: $!";
    # These will count the contigs and the DNA.
    my ($contigs, $bases) = (0, 0);
    # These are used to build the FASTA data.
    my ($contigID, $comment, @rows, @triples);
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(fastaLineIn => 1);
        if ($line =~ /^>(\S+)(?:\s+(.+))?/) {
            # Here we are starting a new contig. Flush the old one and see if we want to keep the
            # new one.
            $stats->Add(fastaContig => 1);
            my ($newID, $newComment) = ($1, $2);
            if ($contigID) {
                push @triples, [$contigID, $comment, join("", @rows)];
            }
            $contigID = $newID;
            $comment = $newComment || '';
            @rows = ();
            $contigs++;
        } else {
            # Here we have a line of DNA data in a contig we're keeping.
            push @rows, $line;
            chomp $line;
            $bases += length($line);
        }
    }
    close $ih;
    # Compute the name.
    my ($name) = $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(id) = ?', [$taxonID], 'scientific-name');
    if (! $name) {
        die "Taxonomic grouping $taxonID not found.";
    } else {
        $name = "$name (modified)";
    }
    my %options = ('sleep' => $opt->sleep, user => $opt->user, password => $opt->password);
    my $gto = RASTlib::Annotate(\@triples, $taxonID, $name, %options);
    $stats->Add(gtoCreated => 1);
    # Output the GTO.
    my $newID = $gto->{id};
    print "Creating new GenomePackage $newID.\n";
    my $newDir = "$outDir/$newID";
    if (-d $newDir) {
        die "Redundant genome ID returned $newID.";
    } else {
        File::Copy::Recursive::pathmk($newDir) || die "Could not create $newDir: $!";
        File::Copy::Recursive::fcopy("$inDir/$file", "$newDir/bin.fa") || die "Could not copy FASTA: $!";
        SeedUtils::write_encoded_object($gto, "$newDir/bin.gto");
        open(my $oh, '>', "$newDir/data.tbl") || die "Could not open new data.tbl: $!";
        print $oh "Genome Name\t$name\n";
        print $oh "Source Package\t$genomeID\n";
        print $oh "Contigs\t$contigs\n";
        print $oh "Base pairs\t$bases\n";
    }
}
print "All done.\n" . $stats->Show();