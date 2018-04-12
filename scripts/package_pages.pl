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
use GPUtils;
use BinningReports;

=head1 Generate Report Web Pages for Genome Packages

    package_pages.pl [ options ] pDir

This script will run through all of the genome packages in a directory. For those with completed quality
evaluations, it will create an HTML page (report.html) that contains a quality report. This page can be accessed to
analyze the problematic roles on PATRIC.

=head2 Parameters

The positional parameter is the name of the genome package directory.

The command-line options are the following.

=over 4

=item roleFile

The C<roles.in.subsystems> file containing the subsystem roles. This is a tab-delimited file: for each role, it contains
a single line with (0) the role ID, (1) the role checksum, and (2) the role name. The default is to use the global SEEDtk
file.

=item missing

If specified, reports will not be generated for packages that already have them.

=item template

The file name of the web page template. The default is C<BinningReports/details.tt> in the kernel library directory.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('pDir',
        ['roleFile|rolefile|r=s', 'role mapping file', { default => "$FIG_Config::global/roles.in.subsystems" }],
        ['missing', 'only process packages without pre-existing reports'],
        ['template=s', 'detail template file', { default => "$FIG_Config::mod_base/kernel/lib/BinningReports/details.tt" }],
        );
my $stats = Stats->new();
# Get the package directory.
my ($pDir) = @ARGV;
if (! $pDir) {
    die "No package directory specified.";
} elsif (! -d $pDir) {
    die "Package directory $pDir missing or invalid.";
}
# Get the missing option.
my $missing = $opt->missing;
# Create the role maps.
my (%cMap, %nameMap);
open(my $rh, '<', $opt->rolefile) || die "Could not open role file: $!";
while (! eof $rh) {
    my ($id, $cksum, $name) = ScriptUtils::get_line($rh);
    $cMap{$cksum} = $id;
    $nameMap{$id} = $name;
    $stats->Add(roleIn => 1);
}
# Load the template.
open(my $ih, '<', $opt->template) || die "Could not open template file: $!";
my $detailsT = join("", <$ih>);
close $ih;
# Build the HTML prefix.
my $prefix = <<'END_HTML';
<html>
<head>
<style type="text/css">
    table.p3basic,
        th, td, tr {
            border-style: double;
            border-collapse: collapse;
            vertical-align: top;
}
</style>
</head>
<body>

END_HTML
my $suffix = "\n</body></html>\n";
# Get all of the packages.
print "Reading packages from $pDir.\n";
my $genomeHash = GPUtils::get_all($pDir);
my $total = scalar(keys %$genomeHash);
print "$total packages found in $pDir.\n";
my $count = 0;
# Loop through the packages.
for my $genome (sort keys %$genomeHash) {
    my $gDir = $genomeHash->{$genome};
    $count++;
    my $skip;
    if (! -s "$gDir/EvalByCheckG/evaluate.log") {
        $stats->Add(packageNoCheckG => 1);
        $skip = 1;
    }
    if (! -s "$gDir/EvalBySciKit/evaluate.log") {
        $stats->Add(packageNoSciKit => 1);
        $skip = 1;
    }
    if ($missing && -s "$gDir/report.html") {
        $stats->Add(skippedAlreadyDone => 1);
        $skip = 1;
    }
    if (! $skip) {
        # Here we want to create the report.
        my $gto = GPUtils::gto_of($genomeHash, $genome);
        $stats->Add(gtoRead => 1);
        print "Analyzing $genome in $gDir ($count of $total).\n";
        BinningReports::UpdateGTO($gto, "$gDir/EvalBySciKit", "$gDir/EvalByCheckG", \%cMap);
        print "Producing HTML.\n";
        my $html = BinningReports::Detail({}, undef, \$detailsT, $gto, \%nameMap);
        open(my $oh, ">$gDir/report.html");
        print $oh "$prefix\n$html\n$suffix\n";
        $stats->Add(reportOut => 1);
        close $oh;
    }
}
print "All done.\n" . $stats->Show();