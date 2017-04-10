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

    parse_bin_report.pl [ options ] inDirectory

This script parses bins.report.txt files and adds the species information to the input file.
The input file should have a sample ID in the first column and a bin number in the second. The bins.report.txt
files must be in subdirectories of the input directory. Each subdirectory should have the same name as a binning
sample.

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
    my @samples = grep { -f "$inDirectory/$_/bins.report.txt" } readdir $dh;
    if (! @samples) {
        die "No completed samples found in $inDirectory.";
    } else {
        # The species datas will go in here. It is a two-level hash, the first key being sample ID and the
        # second the bin number. Each datum will be a list [species, count].
        my %species;
        # Loop through the samples.
        for my $sampleID (@samples) {
            open(my $fh, "<$inDirectory/$sampleID/bins.report.txt") || die "Could not open $sampleID report file: $!";
            while (! eof $fh) {
                my $line = <$fh>;
                if ($line =~ /^BIN\s+(\d+)/) {
                    my $bin = $1;
                    # Loop through the species list. We use an if-undefined operator to keep the
                    # first one.
                    my ($species1, $count) = (undef, 0);
                    while (substr($line, 0, 3) ne '***') {
                        if ($line =~ /^\s+\d+\.\d+:\s+(.+)/) {
                            $species1 //= $1;
                            $count++;
                        }
                        $line = <$fh>;
                    }
                    $species{$sampleID}{$bin} = [$species1, $count];
                }
            }
        }
        # Open the input file.
        my $ih = ScriptUtils::IH($opt->input);
        # Loop through the input.
        while (! eof $ih) {
            my $line = <$ih>;
            my ($sampleID, $bin, @fields) = split /\t/, $line;
            if (@fields) {
                $fields[$#fields] =~ s/[\r\n]+$//;
            }
            my $data = $species{$sampleID}{$bin};
            if ($data) {
                print join("\t", $sampleID, $bin, @fields, @$data) . "\n";
            }
        }
    }
}