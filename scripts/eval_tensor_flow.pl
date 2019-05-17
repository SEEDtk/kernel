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
use File::Copy::Recursive;
use SeedUtils;

=head1 Evaluate a GTO Using Tensor Flow

    eval_tensor_flow.pl [ options ] gto outDir

This script applies the Tensor Flow evaluation tool to a L<GenomeTypeObject>. The output files will be sent to the
specified directory.

=head2 Parameters

There are two positional parameters-- the input GTO file and the output directory. If the output directory exists,
it will be overwritten.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('gto outDir');
# Process the positional parameters.
my ($gto, $outDir) = @ARGV;
if (! $gto) {
    die "No input file specified.";
} elsif (! -s $gto) {
    die "Invalid GTO file $gto.";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir.";
} else {
    File::Copy::Recursive::pathempty($outDir) || die "Error clearing output directory $outDir.";
}
# Create the input file.
print "Creating counts.\n";
my $cmd = "gto_to_roles --counts $gto $FIG_Config::p3data/roles.in.subsystems >$outDir/roles.counts 2>$outDir/roles.not.mapped";
SeedUtils::run($cmd);
# Process it.
print "Computing predictions.\n";
$cmd = "/homes/rbutler/ross/worknn/tfpredict.sh $outDir/roles.counts >$outDir/evaluate.tbl 2>$outDir/err.log";
SeedUtils::run($cmd);
# Evaluate the results.
print "Evaluating results.\n";
$cmd = "process_evaluation --input=$outDir/evaluate.tbl >$outDir/evaluate.log";
SeedUtils::run($cmd);
