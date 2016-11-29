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
use File::Copy::Recursive;
use SeedUtils;

=head1 Run Quality Checks on Packages

    package_check.pl [ options ] dir

This script runs through the packages filling in missing quality reports. Currently, this includes SciKit, TensorFlow,
and CheckM.

=head2 Parameters

The positional parameter is the name of the directory containing the genome packages.

The command-line options are as follows.

=over 4

=item force

All of the evaluations will be performed, even if they already exist.

=item status

If specified, no evaluations will be performed, only the status will be displayed.
=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir package',
        ["force", 'force regeneration of quality data'],
        ["status", 'only show totals'],
        );
my $stats = Stats->new;
# Get the directory and the package.
my ($dir) = @ARGV;
if (! $dir) {
    die "No packages directory specified.";
} elsif (! -d $dir) {
    die "Invalid package directory $dir.";
} else {
    # This is the force flag. We default to the option value.
    my $force = $opt->force;
    # Compute the list of packages to process.
    opendir(my $dh, $dir) || die "Could not open package directory: $!";
    my @packages = sort grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
    # Loop through the packages..
    for my $package (@packages) {
        $stats->Add(packages => 1);
        print "Checking $package.\n";
        my $ok = 1;
        my $pDir = "$dir/$package";
        # Process CheckM.
        my $outDir = "$pDir/EvalByCheckm";
        my $cmd = "checkm lineage_wf --tmpdir $FIG_Config::temp -x fa --file $pDir/evaluate.log $pDir $outDir";
        $ok = Process(CheckM => $outDir, $force, $package, $cmd, $opt->status);
        if ($ok) {
            File::Copy::Recursive::fmove("$pDir/evaluate.log", "$pDir/EvalByCheckm/evaluate.log");
        }
        # Process Tensor Flow.
        $outDir = "$pDir/EvalByTF";
        $cmd = "eval_tensor_flow $pDir/bin.gto $outDir";
        $ok = Process("Tensor Flow" => $outDir, $force, $package, $cmd, $opt->status);
        # Process SciKit.
        $outDir = "$pDir/EvalBySciKit";
        $cmd = "gto_consistency $pDir/bin.gto $outDir $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $FIG_Config::global/roles.to.use";
        $ok = Process("SciKit" => $outDir, $force, $package, $cmd, $opt->status);
    }
}
print "All Done.\n" . $stats->Show();

sub Process {
    my ($type, $outDir, $force, $package, $cmd, $statusOnly) = @_;
    my $retVal = 1;
    if (-d $outDir) {
        if ($statusOnly) {
            $retVal = 0;
            $stats->Add("found-$type" => 1);
        } else {
            if ($force) {
                print "Erasing old $outDir.\n";
                File::Copy::Recursive::pathempty($outDir);
                $stats->Add("clear-$type" => 1);
            } else {
                print "$type already run for $package-- skipping.\n";
                $stats->Add("skip-$type" => 1);
                $retVal = 0;
            }
        }
    }
    if ($retVal) {
        print "Running $type for $package.\n";
        $stats->Add("run-$type" => 1);
        &SeedUtils::run($cmd);
    }
    return $retVal;
}