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

=head1 Describe Script Here

    community_analysis.pl directory

This script analyzes the output from a parameter run on L<community_pipeline.pl> and produces a table of the scores.
For each C<bins.summary.>I<n> file, it will list the parameters used followed by the statistics output. The result will
be in a tab-delimited form that can be loaded into a spreadsheet.

=head2 Parameters

The single positional parameter is the name of the directory containing the community pipeline output files.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory');
# Get the input directory.
my $dir = $ARGV[0];
if (! $dir) {
    die "No directory specified.";
} elsif (! -d $dir) {
    die "Invalid directory $dir.";
} else {
    opendir (my $dh, $dir) || die "Cannot open directory $dir: $!";
    my @summaries = grep { $_ =~ /^bins\.summary\.\d+$/ } readdir ($dh);
    closedir $dh;
    # We now have a list of all the summary files in the directory. We will process each one individually.
    # We build two hashes-- a parm hash and a statistics hash. Each key maps to a subhash of all the values
    # found, keyed by file index.
    my (%parms, %stats, @idxes);
    for my $summaryFile (@summaries) {
        # Compute the corresponding parameter file name.
        my (undef, undef, $idx) = split /\./, $summaryFile;
        push @idxes, $idx;
        my $parmFile = "parms.$idx";
        # Read the parms.
        open (my $ih, "<$dir/$parmFile") || die "Could not open parameter file $idx: $!";
        while (! eof $ih) {
            my $line = <$ih>;
            if ($line =~ /^(\w+)\s+=\s+(.+)/) {
                $parms{$1}{$idx} = $2;
            }
        }
        close $ih;
        undef $ih;
        # Read the statistics from the end of the summary file.
        open ($ih, "<$dir/$summaryFile") || die "Could not open summary file $idx: $1";
        while (! eof $ih && <$ih> !~ /^########/) {};
        while (! eof $ih) {
            my $line = <$ih>;
            if ($line =~ /^(\w+)\s+(\d+)/) {
                $stats{$1}{$idx} = $2;
            }
        }
        close $ih;
    }
    # Now output the table. Note we have a blank column between the parameters and the statistics. This blank column's
    # heading also tells us when to switch hashes.
    my @idxList = sort { $a <=> $b } @idxes;
    my @heads = (sort keys %parms);
    push @heads, '';
    push @heads, sort keys %stats;
    print join("\t", '#', @heads) . "\n";
    for my $idx (@idxList) {
        my @output;
        my $hash = \%parms;
        for my $head (@heads) {
            if (! $head) {
                $hash = \%stats;
                push @output, '';
            } else {
                push @output, ($hash->{$head}{$idx} // '');
            }
        }
        print join("\t", $idx, @output) . "\n";
    }
}
