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

=head1 Print Table of Bins

    bin_table.pl inDirectory

This script reads the FASTA files in a directory and outputs a table of the contig IDs in each file. The table
will be in the form of a tab-delimited file, with the first column being a contig ID and the second being a FASTA
file name.

=head2 Parameters

There is one positional parameter-- the name of the directory containing the input FASTA files.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDirectory');
# Open the input directory and get the FASTA files.
my ($inDirectory) = @ARGV;
if (! $inDirectory) {
    die "Input directory name is required.";
} elsif (! -d $inDirectory) {
    die "Invalid input directory name $inDirectory.";
}
opendir(my $dh, $ARGV[0]) || die "Could not open input directory: $!";
my @files = grep { substr($_,0,1) ne '.' && -f "$inDirectory/$_" } readdir $dh;
closedir $dh;
for my $file (@files) {
    open(my $ih, "<", "$inDirectory/$file") || die "Could not open FASTA file $file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^>(\S+)/) {
            print "$1\t$file\n";
        }
    }
}