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

=head1 Create Taxonomic Profile File for Completeness Computer

    p3x-taxon-profiles.pl [options]

This script is a quick and dirty utility that separates the good genomes in PATRIC by domain.  It is used to create
a repgen-analysis file of last resort that can be used to catch outlier genomes.  The resulting file is fed into
L<p3x-merge-rep-files.pl> with a similarity score of 0.  The standard input should be the C<patric.good.tbl> file,
which contains genome IDs, names, and a full taxonomic lineage for each good genome.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

Additional command-line options are as follows.

=over 4

=item verbose

Display progress messages on STDERR.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::ih_options(),
    ['verbose|debug|v', 'write progress messages to STDERR']);
my $debug = $opt->verbose;
# Set up the counters.
my $stats = Stats->new();
my $gCount = 0;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Skip the header line and write the output header.
my $line = <$ih>;
P3Utils::print_cols(['genome_id', 'name', 'tax_id', 'score']);
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    unless ($line =~ /(\S+)\t([^\t]+)\t\d+::(\d+)/) {
    	print STDERR "WARNING: invalid input line $line" if $debug;
    	$stats->Add(badLine => 1);
    } else {
    	$stats->Add(genomeIn => 1);
    	my ($genome, $name, $tax_id) = ($1, $2, $3);
    	$stats->Add("genome-$tax_id" => 1);
    	$gCount++;
    	P3Utils::print_cols([$genome, $name, $tax_id, 100]);
    	print STDERR "$gCount genomes processed.\n" if $debug && $gCount % 5000 == 0;
    }
}
print STDERR "All done.\n" . $stats->Show() if $debug;