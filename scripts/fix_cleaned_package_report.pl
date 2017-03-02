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

=head1 Modify a Package Report With Trimmed Bins

    fix_cleaned_package_report.pl [ options ] packageReport <trimlogfile >newreport

This script modifies a package report (produced by L<package_report.pl>) to make it possible to relate a trimmed bin
to its original bin. The trimmed bin is given the sample ID from the original bin and its reference genome name is
set to the original bin's name with a modifier of C<coarse> or C<fine>.

The basic procedure is to read in the package report and then run through the
log file (which comes in via the standard input). The information in the log file enables us to connect each
trimmed bin to its original.

=head2 Parameters

The positional parameter is the original package report.

The command-line options are those found in L<ScriptUtils/ih_options> (to specify the incoming log file) plus the following.

=over 4

=item coarse

If specified, it is presumed coarse trimming was used. The default is to presume fine trimming.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageReport', ScriptUtils::ih_options(),
        ['coarse', 'coarse trimming used'],
        );
# Verify the parameter.
my ($packageReport) = @ARGV;
if (! $packageReport) {
    die "No package report specified.";
} elsif (! -s $packageReport) {
    die "Invalid or missing package report file $packageReport.";
}
# Read the package report into a hash.
my %report;
open(my $rh, '<', $packageReport) || die "Could not open package report: $!";
while (! eof $rh) {
    my $line = <$rh>;
    chomp $line;
    my ($sample, $binID, @data) = split /\t/, $line;
    $report{$binID} = [$sample, $binID, @data];
}
close $rh;
# Determine the report type.
my $type = ($opt->coarse ? 'coarse' : 'fine');
# This will be the source bin ID for the current section of the log.
my $sourceBin;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^Processing (\d+\.\d+)/) {
        $sourceBin = $1;
    } elsif ($line =~ /^Creating new GenomePackage (\d+\.\d+)/) {
        # Here we have a new bin created from another one. Fix the sample ID and the reference genome name.
        my $newBin = $1;
        my $sourceLine = $report{$sourceBin};
        my $newLine = $report{$newBin};
        $newLine->[0] = $sourceLine->[0];
        $newLine->[10] = $sourceLine->[2] . " ($type-trimmed)";
    }
}
# Write the report back out.
for my $bin (sort keys %report) {
    my $line = $report{$bin};
    print join("\t", @$line) . "\n";
}
