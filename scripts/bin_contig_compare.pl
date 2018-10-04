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
use Bin;
use File::Copy::Recursive;
use Stats;

=head1 Compare Contig Populations of Bins

    bin_contig_compare.pl [ options ] dir1 dir2 outDir

This script compares the binning results in two sample directories against which a binning job has been run. It is
assumed the same contigs are used in both samples. For each pair of bins (one from the first directory and one from
the second), a list of the contigs in common, the contigs only in the first, and the contigs only in the second will
be produced. A summary report will count the contigs matching and different.

=head2 Parameters

The positional parameters are the names of the two input directories and the name of the output directory.
The output directory will contain a detail report for each pair of bins listing the contigs only in the first,
only in the second, and in both, and a summary report with the overall counts. The summary report will be in
the C<summary.tbl> file and the detail report in a file named I<bin1>C<.>I<bin2>C<.tbl> file, using the
taxon IDs of the bins.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir1 dir2 outDir',
        );
# Verify the parameters.
my ($dir1, $dir2, $outDir) = @ARGV;
if (! $dir1) {
    die "No first input directory specified.";
} elsif (! -d $dir1) {
    die "Invalid or missing directory $dir1.";
} elsif (! -s "$dir1/bins.rast.json") {
    die "$dir1 does not appear to contain a binning job.";
}
if (! $dir2) {
    die "No second input directory specified.";
} elsif (! -d $dir2) {
    die "Invalid or missing directory $dir2.";
} elsif (! -s "$dir2/bins.rast.json") {
    die "$dir2 does not appear to contain a binning job.";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (-f $outDir) {
    die "Invalid output directory $outDir.";
} elsif (! -d $outDir) {
    print "Creaating directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir);
} else {
    print "Output will be to $outDir.\n";
}
# Create the statistics object.
my $stats = Stats->new();
# Open the two binning jobs.
print "Reading bins for $dir1.\n";
my $bins1 = Bin::ReadBins("$dir1/bins.rast.json");
print "Reading bins for $dir2.\n";
my $bins2 = Bin::ReadBins("$dir2/bins.rast.json");
# The summary report will be written here.
open(my $oh, ">$outDir/summary.tbl") || die "Could not open summary.tbl: $!";
print $oh join("\t", qw(Bin1 Bin2 1only both 2only)) ."\n";
# Loop through the bin pairs. We first generate the three lists of contigs, then write the detail report.
for my $bin1 (@$bins1) {
    my $bin1Name = $bin1->taxonID;
    my @bin1Contigs= sort $bin1->contigs;
    my %bin1Contigs = map { $_ => 1 } @bin1Contigs;
    for my $bin2 (@$bins2) {
        my $bin2Name = $bin2->taxonID;
        print "Comparing $bin1Name to $bin2Name.\n";
        # The three contig lists go in here.
        my (@only1, @only2, %both);
        # Find out what occurs in both or only in bin 2.
        for my $contig2 (sort $bin2->contigs) {
            if ($bin1Contigs{$contig2}) {
                $both{$contig2} = 1;
                $stats->Add(both => 1);
            } else {
                push @only2, $contig2;
                $stats->Add(only2 => 1);
            }
        }
        # Find out what occurs only in bin 1.
        for my $contig1 (@bin1Contigs) {
            if (! $both{$contig1}) {
                push @only1, $contig1;
                $stats->Add(only1 => 1);
            }
        }
        # Produce the summary line.
        my @both = sort keys %both;
        my $only2K = scalar @only2;
        my $only1K = scalar @only1;
        my $bothK = scalar @both;
        print $oh join("\t", $bin1Name, $bin2Name, $only1K, $bothK, $only2K) . "\n";
        # Produce the detail report.
        open(my $dh, ">$outDir/$bin1Name.$bin2Name.tbl") || die "Could not open detail file: $!";
        print $dh join("\t", "$dir1.$bin1Name", 'Both', "$dir2.$bin2Name") . "\n";
        while (@both || @only1 || @only2) {
            my $only1 = (shift @only1) // '';
            my $both = (shift @both) // '';
            my $only2 = (shift @only2) // '';
            print $dh join("\t", $only1, $both, $only2) . "\n";
        }
    }
}
print "All done:\n" . $stats->Show();