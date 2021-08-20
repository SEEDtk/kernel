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

=head1 Compute Operon Matrix

    operon_matrix.pl [ options ] inFile outFile

This program computes the operon matrix from a Threonine difference table.  For each row of the table, every pair of
operon values is examined.  If both values are nonempty and nonzero, the ratio will be computed.  At the end, the
mean ratio of each pair will be output.

=head2 Parameters

The positional parameters are the name of the input file and the name of the output file.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inFile outFile',
        );
# Create a statistics object.
my $stats = Stats->new();
my ($inFile, $outFile) = @ARGV;
print "Opening input file $inFile.\n";
open(my $ih, '<', $inFile) || die "Could not opne input file: $!";
my $line = <$ih>;
chomp $line;
my (undef, @operons) = split /\t/, $line;
# We will accumulate results in here.
my %totals;
my %counts;
print "Looping through input file.\n";
while (! eof $ih) {
    $line = <$ih>;
    chomp $line;
    $stats->Add(lineIn => 1);
    my (undef, @prodValues) = split /\t/, $line;
    for (my $i = 0; $i < @prodValues; $i++) {
        my $prodI = $prodValues[$i];
        if ($prodI && $prodI > 0.0) {
            for (my $j = $i + 1; $j < @prodValues; $j++) {
                my $prodJ = $prodValues[$j];
                if ($prodJ && $prodJ > 0.0) {
                    $stats->Add(pairFound => 1);
                    my $ratio = $prodI / $prodJ;
                    my $pair = "$operons[$i]/$operons[$j]";
                    $totals{$pair} += $ratio;
                    $counts{$pair}++;
                }
            }
        }
    }
}
close $ih;
# Write the means to the output.
print "Computing means.\n";
open(my $oh, '>', $outFile) || die "Could not open output file: $!";
print $oh "pair\tmean\tcount\n";
for my $pair (sort keys %totals) {
    $stats->Add(pairOut => 1);
    print $oh join("\t", $pair, $totals{$pair}/$counts{$pair}, $counts{$pair}) . "\n";
}
close $oh;
print "All done.\n" . $stats->Show();
