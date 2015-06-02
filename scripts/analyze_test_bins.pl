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

=head1 Analyze Composition of Community Pipeline Bins Made from Test Contigs

    analyze_bins.pl [ options ] > output.txt

This script reads a bins file and lists the source genomes found in a bin. The input contigs for the bins
file must have been produced by the L<generate_sample.pl> script. Such contigs have an ID that begins with
the ID of the source genome. For each bin, we will list the number of contigs in the bin for each source
genome in it.

=head2 Parameters

The command-line options are those found in L<ScriptUtils/ih_options>, which specifies the input bins file.
The bins file itself is in the format described in L<community_pipeline.pl>, which produces it.

=head2 Output

The output is tab-delimited to the standard output. Each record consists of a genome ID and a contig count.
The data for the bins will be separated by a C<//> indicator.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('parms', ScriptUtils::ih_options(),
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# This hash contains the data for the current bin.
my %binGenomes;
# Loop through the input.
while (! eof $ih) {
    # Get the current line.
    my $line = <$ih>;
    if ($line =~ /^(\d+\.\d+)_/) {
        # Here we have a contig header.
        $binGenomes{$1}++;
    } elsif ($line =~ m#^//#) {
        # Here we have an end-of-bin.
        for my $genome (keys %binGenomes) {
            print "$genome\t$binGenomes{$genome}\n";
        }
        print "//\n\n";
        %binGenomes = ();
    }
}
