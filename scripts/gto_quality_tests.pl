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

=head1 Produce Quality Report on GTOs in a Directory

    gto_quality_tests.pl [ options ] gtoDir

This script will run through all the GTO files in a directory and run the EvalG and EvalCon quality tests on them,
producing a report.

The report will be output to the file C<report.tbl> in the input directory. Status messages will be written to
the standard output.

=head2 Parameters

The sole positional parameter is the name of the input directory containing the GTO files.

The command-line options are the following.

=over 4

=item predictors

The directory containing the function predictors. If omitted, the default function predictors are used.

=item tempDir

The name of a temporary directory to contain working and output files. This directory will be cleared during processing.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('gtoDir',
        ["predictors|P=s", 'function predictors directory for evalCon'],
        ["tempDir|t=s", "working directory for temporary files", { default => "$FIG_Config::temp/QTest" }],
        );
# Get the input directory.
my ($gtoDir) = @ARGV;
if (! $gtoDir) {
    die "No input directory specified.";
} elsif (! -d $gtoDir) {
    die "Input directory $gtoDir not found or invalid.";
}
# Check the function predictors.
my ($predictors, $roleFile, $rolesToUse) = ("$FIG_Config::global/FunctionPredictors", "$FIG_Config::global/roles.in.subsystems",
        "$FIG_Config::global/roles.to.use");
my $funPred = $opt->predictors;
if ($funPred) {
    $predictors = $funPred;
    if (-s "$predictors/roles.in.subsystems") {
        $roleFile = "$predictors/roles.in.subsystems";
    }
    if (-s "$predictors/roles.to.use") {
        $rolesToUse = "$predictors/roles.to.use";
    }
}
print "Function predictors are $predictors with role file $roleFile using $rolesToUse.\n";
# Get the temp directory.
my $tempDir = $opt->tempdir;
if (! -d $tempDir) {
    print "Creating $tempDir.\n";
    File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
} else {
    print "Clearing $tempDir.\n";
    File::Copy::Recursive::pathempty($tempDir) || die "Could not clear $tempDir: $!";
}
# Create the output file.
open(my $oh, ">$gtoDir/report.tbl") || die "Could not open report.tbl: $!";
print $oh join("\t", qw(fileName completeness contamination multiplicity group coarse fine)) . "\n";
# Get the list of GTOs.
print "Reading input directory.\n";
opendir(my $dh, $gtoDir) || die "Could not open $gtoDir: $!";
my @gtos = grep { $_ =~ /\.gto$/ } readdir $dh;
closedir $dh;
print scalar(@gtos) . " found in $gtoDir.\n";
for my $gto (@gtos) {
    # This will be the output line.
    my @cols = ($gto);
    # Get the GTO and run EvalG.
    my $gtoName = "$gtoDir/$gto";
    my $cmd = "check_gto --eval --quiet $tempDir $gtoName";
    print "Executing $cmd.\n";
    SeedUtils::run($cmd);
    # We keep all of the output values from evalG, but we skip over the header line.
    open(my $ih, "<$tempDir/evaluate.log") || die "Could not open EvalG output: $!";
    my $line = <$ih>;
    push @cols, ScriptUtils::get_line($ih);
    close $ih; undef $ih;
    # Run EvalCon.
    $cmd = "gto_consistency $gtoName $tempDir/SciKit $predictors $roleFile $rolesToUse";
    print "Executing $cmd.\n";
    SeedUtils::run($cmd);
    open($ih, "<$tempDir/SciKit/evaluate.log") || die "Could not open EvalCon output: $!";
    my ($fine, $coarse) = (0, 0);
    while (! eof $ih) {
        $line = <$ih>;
        if ($line =~ /^Fine_Consistency=\s+(.+)\%/) {
            $fine = $1;
        } elsif ($line =~ /^Coarse_Consistency=\s+(.+)\%/) {
            $coarse = $1;
        }
    }
    push @cols, $coarse, $fine;
    close $ih; undef $ih;
    print $oh join("\t", @cols) . "\n";
    # Clear the temp directory.
    File::Copy::Recursive::pathempty($tempDir) || die "Could not clear $tempDir: $!";
}
print "All done.\n";
