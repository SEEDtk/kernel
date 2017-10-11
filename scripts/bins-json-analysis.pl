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
use SeedUtils;

=head1 Analyze bins.json File

    bins-json-analysis.pl [ options ] fileName

This script reads a C<bins.rast.json> or C<bins.json> file from a binning run and writes out the ID, name, length, contigs, and reference genome list for each bin.
It is a quick way to find which bins are the result of genome merges.

=head2 Parameters

The positional parameter is the name of the json file.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('fileName',
        );
# Open the input file.
my ($fileName) = @ARGV;
if (! $fileName) {
    die "No file name specified.";
} elsif (! -f $fileName) {
    die "File $fileName not found or invalid.";
}
my $binList = SeedUtils::read_encoded_object($fileName);
for my $bin (@$binList) {
    my $name = $bin->{name};
    my $taxonID = $bin->{taxonID};
    my $len = $bin->{len};
    my $count = scalar @{$bin->{contigs}};
    my $refs = join(", ", @{$bin->{refGenomes}});
    print join("\t", $taxonID, $name, $count, $len, $refs) . "\n";
}
