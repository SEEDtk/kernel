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

=head1 Compare Two Binning Directories

    compare_bins.pl [options] binDir1 binDir2

This script compares the binning results in two directories.  For each sample in both directories,
it will display the number of good and bad bins in both versions.

=head2 Parameters

The positional parameters are the names of the two binning directories.

=cut

use strict;
use P3DataAPI;
use P3Utils;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('binDir1 binDir2');
# Analyze the input directories.
my ($binDir1, $binDir2) = @ARGV;
if (! $binDir2) {
    die "Must specify two directories.";
} elsif (! -d $binDir1) {
    die "Directory $binDir1 missing or invalid.";
} elsif (! -d $binDir2) {
    die "Directory $binDir2 missing or invalid.";
}
# Get all the samples in the first directory.
my $samples1H = getSamples($binDir1);
my $samples2H = getSamples($binDir2);
# Write the comparison.
P3Utils::print_cols(['sample', "$binDir1 good", "$binDir1 bad", "$binDir2 good", "$binDir2 bad"]);
for my $sample (sort keys %$samples1H) {
    my $file2 = $samples2H->{$sample};
    if ($file2) {
        my ($good1, $bad1) = countBins($samples1H->{$sample});
        my ($good2, $bad2) = countBins($file2);
        P3Utils::print_cols([$sample, $good1, $bad1, $good2, $bad2]);
    }
}

## Count the good and bad bins in an index file.
sub countBins {
    my ($file) = @_;
    open(my $ih, '<', $file) || die "Could not open $file: $!";
    # Skip the header line.
    my $line = <$ih>;
    my ($good, $bad) = (0, 0);
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /1$/) {
            $good++;
        } else {
            $bad++;
        }
    }
    return ($good, $bad);
}

## Find all the completed samples in the specified directory.  Returns a hash that maps sample ID to
## eval file.
sub getSamples {
    my ($binDir) = @_;
    opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
    my @samples = grep { -s "$binDir/$_/Eval/index.tbl" } readdir $dh;
    my %retVal = map { $_ => "$binDir/$_/Eval/index.tbl" } @samples;
    return \%retVal;
}
