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
use Shrub;
use Bin;
use Bin::Analyze;
use File::Copy::Recursive;
use Loader;

=head1 Verify Synthetic Bins

    bins_synthetic_check.pl [ options ] binDir

This program displays a report on the contigs in each bin for a run against synthetic data. For each bin, the source
genome from which the synthetic reads were generated is displayed.

=head2 Parameters

There is one positional parameter-- the directory containing the bin information. The json file containing the bin data structures
must be named C<bins.json>. This is output by L<bins_create.pl>. The original FASTA file for the sample contigs must be in
C<contigs.fasta>. This is output by L<bins_coverage.pl>.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDirectory', Shrub::script_options(),
        );
# Create the loader object and get the statistics.
my $loader = Loader->new();
my $stats = $loader->stats;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Create the analyzer.
my $analyzer = Bin::Analyze->new();
# Read the bins.
my $binList;
my ($binDir) = @ARGV;
if (! $binDir) {
    die "A bins directory is required.";
} elsif (! -d $binDir) {
    die "Bin directory $binDir not found or invalid.";
}
# Compute the bins file name.
my $jsonFileName = "$binDir/bins.json";
print "Reading bins.\n";
$binList = Bin::ReadBins($jsonFileName);
# Read in the contig FASTA file and extract the source genome IDs.
my %contigs2Source;
my $fh = $loader->OpenFasta(sampleContigs => "$binDir/contigs.fasta");
while (my $triple = $loader->GetLine(sampleContigs => $fh)) {
    my ($id, $comment) = @$triple;
    # Extract the source genome ID.
    if ($comment =~ /seed_\d+:(\d+\.\d+):/) {
        my $genomeID = $1;
        $contigs2Source{$id} = $genomeID;
        $stats->Add(contigSourceFound => 1);
    } else {
        print STDERR "WARNING: Source genome not found for $id.\n";
        $stats->Add(contigSourceNotFound => 1);
    }
}
# Loop through the bins.
for (my $id = 0; $id < scalar(@$binList); $id++) {
    my $bin = $binList->[$id];
    $stats->Add(binsProcessed => 1);
    # Tally the genome sources in this bin.
    my %genomesFound;
    for my $contig ($bin->contigs) {
        $stats->Add(contigsProcessed => 1);
        my $genome = $contigs2Source{$contig};
        if (! $genome) {
            $stats->Add(contigNotMatched => 1);
        } else {
            $genomesFound{$genome}++;
        }
    }
    # Display this bin's breakdown.
    $analyzer->BinHeader($bin, \*STDOUT, $id);
    my @genomes = sort { $genomesFound{$b} <=> $genomesFound{$a} } keys %genomesFound;
    for my $genome (@genomes) {
        my ($gname) = $shrub->GetEntityValues(Genome => $genome, 'name');
        print "    $genome\t$genomesFound{$genome}\t$gname\n";
    }
    print "\n";
}
print "All done.\n" . $stats->Show();