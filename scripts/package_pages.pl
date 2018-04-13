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
use Template;
use File::Spec;
use File::Copy::Recursive;

=head1 Generate Report Web Pages for Genome Packages

    package_pages.pl [ options ] pDir1 pDir2 ... pDirN

This script will run through all of the genome packages in one or more directories. For those with completed quality
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

=item tDir

The name of the directory containing the web page templates.

=item wDir

If specified, the name of a web directory to contain a copy of the pages. A master page will be created with the name
C<index.html> in this directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('pDir1 pDir2 ... pDirN',
        ['roleFile|rolefile|r=s', 'role mapping file', { default => "$FIG_Config::global/roles.in.subsystems" }],
        ['missing', 'only process packages without pre-existing reports'],
        ['tDir|tdir|templates=s', 'template file directory', { default => "$FIG_Config::mod_base/kernel/lib/BinningReports" }],
        ['wDir|wdir|webDir|webdir|w=s', 'proposed web directory']
        );
my $stats = Stats->new();
# Get the package directory.
my (@pDirs) = @ARGV;
if (! @pDirs) {
    die "No package directory specified.";
}
# Check for a web directory.
my $wDir = $opt->wdir;
if ($wDir) {
    if (! -d $wDir) {
        print "Creating web directory $wDir.\n";
        File::Copy::Recursive::pathmk($wDir) || die "Could not create $wDir: $!";
    } else {
        print "Clearing web directory $wDir.\n";
        File::Copy::Recursive::pathempty($wDir) || die "Could not empty $wDir: $!";
    }
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
# This will hold the parameters for the master report.
my %master;
# Load the detail template.
my $tDir = $opt->tdir;
open(my $ih, "<$tDir/details.tt") || die "Could not open detail template file: $!";
my $detailsT = join("", <$ih>);
close $ih;
# Build the HTML prefix and suffix.
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
# Loop through the directories.
for my $pDir (@pDirs) {
    if (! -d $pDir) {
        print "WARNING: Directory $pDir not found-- skipping.\n";
    } else {
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
                # Is this a good bin?
                my $details = $gto->{genome_quality_measure};
                my $goodSeed = (GPUtils::good_seed($gto) ? 'Y' : '');
                my $coarse = $details->{consis_data}{Coarse_Consistency} // '';
                my $fine =   $details->{consis_data}{Fine_Consistency} // '';
                my $complt = $details->{checkg_data}{Completeness} // '';
                my $contam = $details->{checkg_data}{Contamination} // '';
                my $good = ($goodSeed && $fine >= 85 && $complt >= 80 && $contam <= 15);
                # Compute the report URL.
                my $reportUrl = "$genome.html";
                # Add our data to the master report.
                my $genomeData = {genome_url => "https://www.patricbrc.org/view/Genome/$genome",
                        genome_id => $genome, report_url => $reportUrl, genome_name => $gto->{scientific_name},
                        scikit_coarse => $coarse, scikit_fine => $fine, checkg_completeness => $complt,
                        checkg_contamination => $contam, good_seed => $goodSeed};
                if ($good) {
                    $master{good_count}++;
                    push @{$master{good}}, $genomeData;
                } else {
                    $master{bad_count}++;
                    push @{$master{bad}}, $genomeData;
                }
                # If there is a web directory, copy the report.
                if ($wDir) {
                    File::Copy::Recursive::fcopy("$gDir/report.html", "$wDir/$reportUrl") || die "Could not copy to $reportUrl: $!";
                }
            }
        }
    }
}
if ($wDir) {
    # Now output the master report.
    my $templateEngine = Template->new(ABSOLUTE => 1);
    my $vars = \%master;
    my $html;
    open(my $th, "<$tDir/master.tt") || die "Could not open master template: $!";
    $templateEngine->process($th, $vars, \$html);
    open(my $oh, ">$wDir/index.html") || die "Could not open master output: $!";
    print $oh "$prefix\n$html\n$suffix\n";
    close $oh;
}
print "All done.\n" . $stats->Show();