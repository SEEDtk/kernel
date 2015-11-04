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
use Bin;


=head1 Compare Community Bins

    bins_compare.pl [ options ] binDir

This program compares a set of binning runs. The contigs are sorted by mean coverage, and then the assigned bin for each run
is listed in columns next to the contig ID.

=head2 Parameters

The positional parameters are the name of the contigs file followed by the json files for the assigned bins.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('contigs bins1 bins2 ...');
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
# Read the bins.
my $binList;
my ($contigFile, @binFiles) = @ARGV;
if (! $contigFile) {
    die "A contig file is required.";
} elsif (! -s $contigFile) {
    die "Contig file $contigFile not found or invalid.";
}
# Read in the contigs.
print STDERR "Processing contig file $contigFile.\n";
my $contigList = Bin::ReadContigs($contigFile);
# Now we process the bins. For each each contig, we are going to keep a list. The list will indicate which bin contains the
# contig in each bin file.
my %contigData = map { $_->contig1 => [] } @$contigList;
$stats->Add(contigs => scalar(@$contigList));
for (my $iBin = 0; $iBin < @binFiles; $iBin++) {
    my $binFile = $binFiles[$iBin];
    print STDERR "Processing bin file $binFile.\n";
    $stats->Add(binFiles => 1);
    my $binList = Bin::ReadBins($binFile);
    for my $bin (@$binList) {
        $stats->Add(binsProcessed => 1);
        # Get the name of the bin.
        my $bName = $bin->contig1;
        # Get its contigs.
        my @contigsInBin = $bin->contigs;
        # Is it a singleton bin?
        if (scalar(@contigsInBin) <= 1) {
            # Yes. Blank the name.
            $bName = 'X';
        }
        # Loop through the contigs in this bin.
        for my $contig ($bin->contigs) {
            $stats->Add(contigsProcessed => 1);
            $contigData{$contig}[$iBin] = $bName;
        }
    }
    # Fill in the blanks.
    for my $contig (keys %contigData) {
        $contigData{$contig}[$iBin] //= ' ';
    }
}
# Write the report.
my @sorted = sort { $b->meanCoverage <=> $a->meanCoverage } @$contigList;
for my $contigBin (@sorted) {
    my $contigID = $contigBin->contig1;
    my $locList = $contigData{$contigID};
    my %binCounts;
    for my $loc (@$locList) {
        $binCounts{$loc}++;
    }
    my ($bestCount) = sort { $b <=> $a } values %binCounts;
    $bestCount //= 0;
    print join("\t", $contigID, $contigBin->meanCoverage, @$locList, $bestCount) . "\n";
}
print STDERR "All done.\n" . $stats->Show();