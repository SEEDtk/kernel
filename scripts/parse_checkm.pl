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


=head1 CheckM Output Formatting

    parse_checkm.pl [ options ] inDirectory

This script parses checkm output files and adds the completeness and contamination scores to the input file.
The input file should have a sample ID in the first column and a bin number in the second. The checkm output
files must be in the form I<sampleID>C<.run.log> in the input directory. The bins should be named C<bin>I<X>,
where I<X> is the bin number.

=head2 Parameters

The single positional parameter is the name of the input directory.

The command-line options are those found in L<ScriptUtils/ih_options>, which specify the standard input.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDirectory', ScriptUtils::ih_options(),
        );
# Process the input directory.
my ($inDirectory) = @ARGV;
if (! $inDirectory) {
    die "No input directory specified.";
} elsif (! -d $inDirectory) {
    die "Invalid input directory specified.";
} else {
    opendir(my $dh, $inDirectory) || die "Could not open $inDirectory: $!";
    my @files = grep { $_ =~ /\.run\.log$/ } readdir $dh;
    if (! @files) {
        die "No checkm run.log files found in $inDirectory.";
    } else {
        # The scores will go in here. It is a two-level hash, the first key being sample ID and the
        # second the bin number.
        my %scores;
        # Loop through the files.
        for my $file (@files) {
            my ($sampleID) = ($file =~ /(.+)\.run\.log$/);
            open(my $fh, "<$inDirectory/$file") || die "Could not open $sampleID log file: $!";
            while (! eof $fh) {
                my $line = <$fh>;
                if ($line =~ /^\s+bin(\d+)\s+(\w+)\s\(UID\d+\)\s+\d+\s+\d+\s+\d+/) {
                    my $bin = $1;
                    my $cat = $2;
                    my @fields = split /\s+/, $line;
                    my $het = pop @fields;
                    my $contam = pop @fields;
                    my $score = pop @fields;
                    my @scoreList = ($score, $contam, $cat);
                    $scores{$sampleID}{$bin} = \@scoreList;
                }
            }
        }
        # Open the input file.
        my $ih = ScriptUtils::IH($opt->input);
        # Loop through the input.
        while (! eof $ih) {
            my $line = <$ih>;
            $line =~ s/[\r\n]+$//;
            my ($sampleID, $bin, @fields) = split /\t/, $line;
            my $scoreList = $scores{$sampleID}{$bin};
            if ($scoreList) {
                print join("\t", $sampleID, $bin, @fields, @$scoreList) . "\n";
            }
        }
    }
}