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

    package_check.pl [ options ] dir pkg1 pkg2 ... pkgN

This script runs through the packages filling in missing quality reports. Currently, this includes SciKit (EvalCon) and
and CheckG (EvalG). CheckM is optionally available.

=head2 Parameters

The first positional parameter is the name of the directory containing the genome packages. Any additional parameters are the
IDs of the genome packages to check. If no IDs are specified, all packages are checked.

The command-line options are as follows.

=over 4

=item force

All of the evaluations will be performed, even if they already exist. A value of C<SciKit> can be specified to only
rerun the SciKit evaluations, and a value of C<CheckG> can be specified to only rerun the CheckG evaluations.

=item status

If specified, no evaluations will be performed, only the status will be displayed.

=item clean

If specified, empty directories will be deleted.

=item checkm

If specified, CheckM evaluations will be performed.

=item predictors

Name of the function predictors directory. The default is FunctionPredictors in the SEEDtk global data directory.
If this option is specified, the role files in the predictors directory will be used instead of the global role files.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir pkg1 pkg2 ... pkgN',
        ["force:s", 'force regeneration of quality data'],
        ["checkm", 'do CheckM evaluations'],
        ["status", 'only show totals'],
        ["clean", 'delete empty packages'],
        ["predictors=s", 'function predictors directory'],
        );
my $stats = Stats->new;
print "Processing options.\n";
my ($rolesToUse, $roleFile, $predictors);
if ($opt->predictors) {
    $predictors = $opt->predictors;
    $rolesToUse = "$predictors/roles.to.use";
    $roleFile = "$predictors/roles.in.subsystems";
} else {
    $predictors = "$FIG_Config::p3data/FunctionPredictors";
    $rolesToUse = "$FIG_Config::p3data/roles.to.use";
    $roleFile = "$FIG_Config::p3data/roles.in.subsystems";
}
print "Scikit files are $predictors, $rolesToUse, and $roleFile.\n";
# Get the directory and the package.
my ($dir, @packages) = @ARGV;
if (! $dir) {
    die "No packages directory specified.";
} elsif (! -d $dir) {
    die "Invalid package directory $dir.";
} else {
    # This is the force flag. We default to the option value.
    my $force = $opt->force;
    my %force;
    if (defined $force) {
        if ($force eq '') {
            $force{SciKit} = 1;
            $force{CheckG} = 1;
        } else {
            $force{$force} = 1;
        }
    }
    # This is the cleaning flag.
    my $clean = $opt->clean;
    # This is the checkM flag. We only do CheckM if it is TRUE.
    my $checkmFlag = ($opt->checkm ? 1 : 0);
    # Compute the list of packages to process.
    if (! @packages) {
        opendir(my $dh, $dir) || die "Could not open package directory: $!";
        @packages = sort grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
        print scalar(@packages) . " directories found.\n";
    } else {
        print "Packages taken from command line.\n";
    }
    # Loop through the packages..
    for my $package (@packages) {
        $stats->Add(packages => 1);
        if (! $opt->status) {
            print "Checking $package.\n";
        }
        my $ok = 1;
        my $pDir = "$dir/$package";
        if (! -s "$pDir/bin.gto") {
            if ($clean) {
                File::Copy::Recursive::pathrmdir($pDir);
                print "Empty directory $package removed.\n";
            } else {
                print "WARNING: package $package is empty!\n";
            }
            $stats->Add(emptyPackages => 1);
        } else {
            my ($outDir, $cmd);
            # Process CheckG.
            $outDir = "$pDir/EvalByCheckG";
            $cmd = "check_gto --eval --quiet $outDir $pDir/bin.gto";
            $ok = Process(CheckG => $outDir, $force{CheckG}, $package, $cmd, $opt->status);
            # Process SciKit.
            $outDir = "$pDir/EvalBySciKit";
            $cmd = "gto_consistency $pDir/bin.gto $outDir $predictors $roleFile $rolesToUse";
            $ok = Process(SciKit => $outDir, $force{SciKit}, $package, $cmd, $opt->status);
            # Process CheckM if the user wants it.
            $outDir = "$pDir/EvalByCheckm";
            $cmd = "checkm lineage_wf --tmpdir $FIG_Config::temp -x fa --file $pDir/evaluate.log $pDir $outDir";
            if ($checkmFlag) {
                $ok = Process(CheckM => $outDir, 0, $package, $cmd, $opt->status);
                if ($ok) {
                    File::Copy::Recursive::fmove("$pDir/evaluate.log", "$pDir/EvalByCheckm/evaluate.log");
                }
            }
        }
    }
}
print "All Done.\n" . $stats->Show();

sub Process {
    my ($type, $outDir, $force, $package, $cmd, $statusOnly) = @_;
    my $retVal = ! $statusOnly;
    if (-d $outDir) {
        if ($statusOnly) {
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
    } elsif (! $statusOnly) {
        # Create the output directory. If we fail, then we have a race condition.
        $retVal = File::Copy::Recursive::pathmk($outDir);
        if (! $retVal) {
            print "Failed to acquire directory for $type of $package-- skipping.\n";
            $stats->Add("race-$type" => 1);
        }
    }
    if ($retVal) {
        print "Running $type for $package.\n";
        $stats->Add("run-$type" => 1);
        &SeedUtils::run($cmd);
    }
    return $retVal;
}