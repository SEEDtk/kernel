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
use Stats;
use Bin::Score;

=head1 Audit Community Bins

    bins_audit.pl [ options ] contigs json

This program performs an audit of the major bins in a binning run, where a major bin is defined here as one with
multiple contigs. For each such bin, it displays the mean and standard deviation of the coverage for the bin, the
coverage of the bin's seed contig, and a list of the missing universal roles.
=head2 Parameters

The two positional parameters are the name of a contigs file (C<contigs.ref.bins> or C<contigs.bins>) and the name of
a JSON file containing the final bins (C<bins.json> or C<bins.new.json>).

The command-line options are the following.

=over 4

=item unifile

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column,
occurrence counts in the second, and role names in the third.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('contigFile binFile',
                ['unifile=s',      'universal role file', { default => "$FIG_Config::global/uni_roles.tbl" }],
        );
# Create the statistics object.
my $stats = Stats->new();
# Get the file names.
my ($contigFile, $binFile) = @ARGV;
if (! $contigFile) {
    die "A contigs file name is required.";
} elsif (! $binFile) {
    die "A bin file name is required.";
}
# Get the universal roles.
my $uniHash = Bin::Score::ReadUniHash($opt->unifile);
# Read the contigs.
my %contigs = map { $_->contig1 => $_ } @{Bin::ReadContigs($contigFile)};
$stats->Add(contigsRead => scalar keys %contigs);
# Read the bins.
my $binList = Bin::ReadBins($binFile);
$stats->Add(binsRead => scalar @$binList);
# Loop through the bins, producing the report.
for my $bin (@$binList) {
    my @contigList = $bin->contigs;
    if (scalar @contigList <= 1) {
        $stats->Add(singletonBin => 1);
    } else {
        $stats->Add(goodBin => 1);
        my $contig1 = $bin->contig1;
        # Get the seed contig coverage.
        my $coverage1 = $contigs{$contig1}->meanCoverage;
        # Loop through the contigs in the bin, computing the mean and standard deviation of the coverage.
        my ($total, $squareT, $n, $cmin, $cmax) = (0, 0, 0, undef, 0);
        for my $contig ($bin->contigs) {
            $stats->Add(contigInBin => 1);
            my $c = $contigs{$contig}->meanCoverage;
            $cmax = $c if $c > $cmax;
            $cmin = $c if ! defined $cmin || $c < $cmin;
            $total += $c;
            $squareT += $c*$c;
            $n++;
        }
        my $m = $total / $n;
        my $sd = sqrt($squareT / $n - $m * $m);
        # Now compute the missing universal roles.
        my $uniBin = $bin->uniProts;
        my @missing;
        for my $role (keys %$uniHash) {
            if (! $uniBin->{$role}) {
                push @missing, $uniHash->{$role};
                $stats->Add(uniRoleMissing => 1);
            } else {
                $stats->Add(uniRoleFound => 1);
            }
        }
        # Print the report on this bin.
        print "BIN $contig1. Coverage = $m (+/- $sd, $n contigs ranging from $cmin to $cmax) from $coverage1.\n";
        print "    Missing Roles\n";
        print "    -------------\n";
        for my $missing (@missing) {
            print "     $missing\n";
        }
        print "\n\n";
    }
}
print "All done.\n" . $stats->Show();