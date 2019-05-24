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

=head1 Delete Samples from a Binning Directory

    p3-delete-samples.pl [options] binDir

This script takes a list of sample IDs as input and deletes the identified samples from the specified binning directory.
Each sample is a subdirectory of the indicated binning directory.

=head2 Parameters

The positional parameter is the name of the binning directory.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The command-line options in <P3Utils/col_options> can be used to identify the input column containing the sample IDs.

=cut

use strict;
use P3Utils;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('binDir', P3Utils::col_options(), P3Utils::ih_options(),
        );
# Get the input directory name.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "Input directory $binDir missing or invalid.";
}
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Get the sample IDs.
my $samples = P3Utils::get_col($ih, $keyCol);
# Loop through the samples.
for my $sample (@$samples) {
    if (! -d "$binDir/$sample") {
        print "$sample is already deleted.\n";
    } else {
        print "Deleting $sample.\n";
        File::Copy::Recursive::pathempty("$binDir/$sample") || die "Could not delete $sample: $!";
        rmdir "$binDir/$sample";
    }
}
