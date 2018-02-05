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

=head1 SRA Sample ID Scraper

    sra_scrape.pl [ options ]

This is a simple script that extracts SRA sample IDs from an NCBI download file and outputs them as a single-column table.

=head2 Parameters

There are no positional parameters. The file reads the standard input and writes to the standard output.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# This will hold the IDs we find.
my %srs;
# Write the output header.
print "sample\n";
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /SRA:\s+(SRS\d+)/) {
        $srs{$1} = 1;
    }
}
# Write the output.
for my $srsID (sort keys %srs) {
    print "$srsID\n";
}
