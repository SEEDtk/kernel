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

=head1 Parse Predictions from Rick and FangFang Models

    parse_cancer_predictions [ options ] rickFile ffFile

This script writes out predictions in the Rick and FangFang models. Both input files are tab-delimited with headers. The Rick input file has
two header records.

The first column in both files contains a CCLE cell line identifier. The Rick file contains three columns for each drug. The first header is the drug ID in all
three columns. The third column (header C<Logical>) is C<1> for a positive response and C<0> for no response. The FangFang file contains one column for each
drug with the drug ID as the column header. A positive number indicates a positive response and a negative number indicates a negative response.

We will output a tab-delimited file with four columns (0) drug name, (1) cell line name, (2) Rick prediction number, (3) FangFang prediction number.
Only combinations present in both files will be output.

=head2 Parameters

The positional parameters are the names of the two files. The Rick file must be first.

The command-line options are as follows:

=over 4

=item filter

The name of a file containing acceptable cell lines in the first column. Only these cell lines will be included.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('rickFile ffFile',
        ['filter|f=s', 'filter file']
        );
# Get the options.
my $filtered = $opt->filter;
my %filter;
if ($filtered) {
    open(my $fh, "<$filtered") || die "Could not open filter file: $!";
    my $line = <$fh>;
    while (! eof $fh) {
        my ($cl) = split /\t/, <$fh>;
        $filter{$cl} = 1;
    }
}
# This hash will map {cellLine}{drug} to the Rick number.
my %rHash;
# This hash will map each result column in the Rick file to its drug ID.
my %rCols;
# Open the Rick file and parse the headers.
open(my $ih, "<$ARGV[0]") || die "Could not open Rick file: $!";
my @headers = ScriptUtils::get_line($ih);
my @head2 = ScriptUtils::get_line($ih);
for (my $i = 0; $i < @headers; $i++) {
    if ($head2[$i] eq 'Prob') {
        $rCols{$i} = $headers[$i];
    }
}
# Loop through the Rick file.
while (! eof $ih) {
    my @data = ScriptUtils::get_line($ih);
    my $cl = $data[0];
    $cl =~ s/^CCLE\.//;
    if ($cl && (! $filtered || $filter{$cl})) {
        for my $i (keys %rCols) {
            $rHash{$cl}{$rCols{$i}} = $data[$i];
        }
    }
}
close $ih; undef $ih;
# Get a hash of the Rick drugs.
my %dHash = map { $rCols{$_} => $rCols{$_} } keys %rCols;
# Compute their names.
for my $type (qw(CCLE CTRP gCSI GDSC NCI60)) {
    open(my $dh, "<$FIG_Config::data/Cancer/drugs/${type}_drugs") || die "Could not open $type drug file: $!";
    my $header = <$dh>;
    while (! eof $dh) {
        my $line = <$dh>;
        if ($line =~ /^(\S+)\t[^\t]+\t(.+)/) {
            if ($dHash{$1}) {
                $dHash{$1} = $2;
            }
        }
    }
}
# This hash will map {cellLine}{drug} to the FangFang indicator.
my %fHash;
# Read the headers.
open($ih, "<$ARGV[1]") || die "Could not open FangFang file: $!";
@headers = ScriptUtils::get_line($ih);
# Loop through the FangFang file.
while (! eof $ih) {
    my @data = ScriptUtils::get_line($ih);
    my $cl = $data[0];
    $cl =~ s/^CCLE\.//;
    if (! $filtered || $filter{$cl}) {
        for (my $i = 1; $i < @data; $i++) {
            $fHash{$cl}{$headers[$i]} = $data[$i];
        }
    }
}
# Get the drugs in common.
my @common = sort grep { $dHash{$_} } @headers;
# Now write the results.
print join("\t", qw(Drug Cell_line Rick FangFang)) . "\n";
for my $cl (sort keys %rHash) {
    my $rData = $rHash{$cl};
    my $fData = $fHash{$cl} // {};
    for my $drug (@common) {
        my $drugName = $dHash{$drug};
        print join("\t", $drugName, $cl, $rData->{$drug} // '', $fData->{$drug} // '') . "\n";
    }
}
