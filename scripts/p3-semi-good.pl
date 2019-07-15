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

=head1 Process Evaluation File to Produce Semi-Good Genome List

    p3-semi-good.pl [options]

This specialized script reads through the C<patric.sort.tbl> file and output the genomes that would be good if they had
a good seed protein.  (This includes the ones that are genuinely good.)  It also converts the taxon lineage to a taxon
ID.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

=cut

use strict;
use GEO;
use P3Utils;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('parms', P3Utils::ih_options());
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read and echo the header line.
my @headers = P3Utils::get_fields($ih);
$headers[2] = "taxon_id";
$headers[7] = "semi_good";
pop @headers;
P3Utils::print_cols(\@headers);
# Loop through the data lines.
while (! eof $ih) {
    my @data = P3Utils::get_fields($ih);
    # Fix the taxon ID.
    if ($data[2] =~ /::(\d+)$/) {
        $data[2] = $1;
    }
    # Check the goodness and write the line.
    if (GEO::consistX($data[4]) && GEO::completeX($data[5]) && GEO::contamX($data[6])) {
        $data[7] = "1";
        pop @data;
        P3Utils::print_cols(\@data);
    }
}
