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

=head1 Prepare for Cross Validation

    cross_validate.pl [ options ] modelDir

This script will set up alternate training files for cross-validation of a dl4j.run model.  It creates
multiple copies of the training file and creates a vparms.prm file to run them all as a search.

=head2 Parameters

The positional parameter is the name of the model directory.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('modelDir'
        );
my ($modelDir) = @ARGV;
if (! $modelDir) {
    die "No model directory specified.";
} elsif (! -s "$modelDir/parms.prm") {
    die "No parms.prm found in $modelDir.";
}
# Open the input and output parms files.
print "Using model directory $modelDir.\n";
open(my $ih, '<', "$modelDir/parms.prm") || die "Could not open input parms.prm: $!";
open(my $oh, '>', "$modelDir/vparms.prm") || die "Could not open output vparms.prm: $!";
# Read the parms file and echo to vparms.prm.
my $inFile = "$modelDir/training.tbl";
my $testSize = 2000;
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^--input\s+(\S+)/) {
        $inFile = $1;
    } else {
        print $oh $line;
        if ($line =~ /^--testSize\s+(\d+)/) {
            $testSize = $1;
        }
    }
}
close $ih;
# This will form the list of input files we're creating.  We'll write them out later.
my @inputs;
# Now we read in the entire training file.
open(my $th, '<', $inFile) || die "Could not open training file $inFile: $!";
# Read the whole training set.
my @trainSet = <$th>;
close $th;
# Step through the training set in testing-set sizes.
my $fileCount = 1;
my $fileSize = scalar @trainSet;
my $train1 = 1;
while ($train1 < $fileSize) {
    open(my $oth, '>', "$modelDir/training$fileCount.tbl") || die "Could not open output file $fileCount: $!";
    print $oth $trainSet[0];
    my $train0 = $train1;
    # Write the testing set.
    for (my $testOut = 0; $testOut < $testSize && $train1 < $fileSize; $testOut++) {
        print $oth $trainSet[$train1];
        $train1++;
    }
    # Write everything but the testing set.
    for (my $i = 1; $i < $fileSize; $i++) {
        if ($i < $train0 || $i >= $train1) {
            print $oth $trainSet[$i];
        }
    }
    close $oth;
    push @inputs, "$modelDir/training$fileCount.tbl";
    print "Training file $fileCount built.\n";
    $fileCount++;
}
print $oh "--input " . join(", ", @inputs) . "\n";
close $oh;
