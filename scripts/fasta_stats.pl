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
use FastA;
use Statistics::Descriptive;

=head1 Compute FASTA Statistics

    fasta_stats.pl [ options ]

This script reads a fasta file and outputs the number of records, the mean sequence length, and the standard
deviation of the sequence length.

=head2 Parameters

There are no positional parameters or command-line options.  The FASTA file should be the standard input.

=cut

$| = 1;
# Open the input file.
my $fh = FastA->new(\*STDIN);
# Initialize the stats object.
my $stats = Statistics::Descriptive::Sparse->new();
while ($fh->next()) {
    my $len = length $fh->left();
    $stats->add_data($len);
}
print $stats->count() . " sequences. Mean = " . $stats->mean() . " with sdev = " . $stats->standard_deviation() . ".\n";
