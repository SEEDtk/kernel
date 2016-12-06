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

=head1 Show Statistics for Binning Package Contigs

    package_contigs.pl [ options ] outFile

After a bin has been converted into a genome package, this script will analyze the contigs to determine the
coverages. The output file will contain a list of the contigs, sorted from longest to shortest, along with
each contig's length and coverage. The contig ID must be in the standard format produced by the SPAdes
assembler. The standard output will show the mean and standard deviation of the coverage.

=head2 Parameters

There is a single positional parameter-- the name of the output file.

The standard input should be the FASTA file for a bin. It can be specified via the command-line options in L<ScriptUtils/ih_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outFile', ScriptUtils::ih_options(),
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# We need a list of contig info, plus the numbers for computing mean and sdev.
my ($sum, $sqr, @contigs) = (0, 0);
# This will be used to compute the mean and standard deviation.
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^>(\S+)/) {
        my $id = $1;
        if ($id =~ /length_(\d+)_cov_(.+)/) {
            my ($len, $cov) = ($1, $2);
            push @contigs, [$id, $len, $cov];
            $sum += $cov;
            $sqr += $cov*$cov;
        }
    }
}
# Open the output file. If there is none, there is no output.
my ($outFile) = @ARGV;
if (! $outFile) {
    print "No output file. Contigs not printed.\n";
} else {
    open(my $oh, '>', $outFile) || die "Could not open output file: $!";
    my @sorted = sort { $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] || $a->[0] cmp $b->[0] } @contigs;
    for my $sorted (@sorted) {
        print $oh join("\t", @$sorted) . "\n";
    }
}
my $n = scalar @contigs;
if (! $n) {
    print "No contigs found or improper contig ID format.\n";
} else {
    my $mean = $sum / $n;
    my $sdev = sqrt($sqr / $n - $mean * $mean);
    print "Mean = $mean. Sdev = $sdev.\n";
}
