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
use SeedUtils qw();
use BinningReports;
use File::Copy::Recursive;

=head1 Test Binning Reports

    binning_report_test.pl [ options ] testDir

This script is used to test the L<BinningReports> module. it takes as input a directory containing key files, then loads the files
into memory and calls the main L<BinningReports/Process> method.

=head2 Parameters

The sole positional parameter should be the name of the test directory. This directory will contain the output HTML files, and should contain the
following input files.

=over 4

=item bins.json

A JSON file containing the bin definition structure.

=item ppr_report.json

A JSON file containing the potentially-problematic-roles report structure.

=item quality.json

A JSON file containing the quality report structure.

=item summary.tt

An HTML template file for the summary report.

=item details.tt

An HTML template file for each detail report.

=item params.json

A JSON file containing the parameter structure.

=back

The structures encoded into JSON files are all described in L<BinningReports>.

The command-line options are the following.

=over 4

=item group

If specified, the genome group path to use. If a genome group path is specified, it changes the report output.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('testDir',
        ['group=s', 'genome group path']
        );
# Use a dummy job ID.
my $jobID = "58ce6a77-e131-40fe-982b-6a8359e0d4e7";
# Get the test directory.
my ($testDir) = @ARGV;
if (! $testDir) {
    die "No test directory specified.";
} elsif (! -d $testDir) {
    die "Test directory invalid or missing.";
} elsif (! -f "$testDir/summary.tt") {
    die "Summary template file not found in $testDir.";
} elsif (! -f "$testDir/details.tt") {
    die "Detail template file not found in $testDir.";
}
# Create the role map.
my %roleMap;
open(my $ih, '<', "$FIG_Config::global/roles.in.subsystems") || die "Could not open role map: $!";
while (my $line = <$ih>) {
    if ($line =~ /^(\S+)\t\S+\t(.+)/) {
        $roleMap{$1} = $2;
    }
}
close $ih; undef $ih;
# Get the input files.
print "Reading input files.\n";
my $bins_json = SeedUtils::read_encoded_object("$testDir/bins.json");
my $ppr_report_json = SeedUtils::read_encoded_object("$testDir/ppr_report.json");
my $quality_json = SeedUtils::read_encoded_object("$testDir/quality.json");
my $params = SeedUtils::read_encoded_object("$testDir/params.json");
# Clear the old reports.
opendir(my $dh, $testDir) || die "Could not open $testDir: $!";
my @old = grep { $_ =~ /^\d+\.\d+$/ && -d "$testDir/$_" } readdir $dh;
close $dh; undef $dh;
for my $dir (@old) {
    print "Deleting $testDir/$dir\n";
    File::Copy::Recursive::pathrmdir("$testDir/$dir") || die "Could not delete $dir: $!";
}
# Generate the reports.
print "Invoking binning reports module.\n";
my ($summary, $detailsH) = BinningReports::Process($jobID, $params, $quality_json, $ppr_report_json, $bins_json, "$testDir/summary.tt", "$testDir/details.tt", $opt->group,
        \%roleMap);
# Write out the reports.
print "Storing output files.\n";
write_html($summary, "$testDir/summary.html");
for my $dir (keys %$detailsH) {
    File::Copy::Recursive::pathmk("$testDir/$dir") || die "Could not create $dir: $!";
    write_html($detailsH->{$dir}, "$testDir/$dir/report.html");
}

sub write_html {
    my ($text, $fileName) = @_;
    print "Writing report to $fileName.\n";
    open(my $oh, ">$fileName") || die "Could not open $fileName: $!";
    print $oh "<head><style>
                table.p3basic {
                    border-collapse: separate;
                    empty-cells: hide;
                    border-bottom: #d3d3d3 1px solid;
                    border-right: #d3d3d3 1px solid;
                    border-spacing: 0;
                }
                table.p3basic th, table.p3basic td, table.p3basic caption {
                    padding: 4px 5px;
                }
                table.p3basic thead th {
                    background: #fff;
                    color: #666666;
                    font-weight: bold;
                }
                table.p3basic td, th {
                    border-top: #d3d3d3 1px solid;
                    border-left: #d3d3d3 1px solid;
                }
               </style><body>\n";
    print $oh $text;
    print $oh "</body>\n";
    close $oh;
}
