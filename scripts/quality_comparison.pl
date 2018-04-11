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
use Stats;

=head1 Compare Quality Files

    quality_comparison.pl [ options ] file1 file2

This script compares the scores in two quality files. Only rows for which all four scores are present in both files will
be compared. The output file contains the bin ID, both sets of scores, and both good-genome flags.

=head2 Parameters

The positional parameters are the names of the two quality files. The first is considered I<old> and the second I<new>.

The output file will be tab-delimited with headers.

Progress messages are written to the standard error output.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('file1 file2',
        );
my $stats = Stats->new();
# Get the file names.
my ($file1, $file2) = @ARGV;
if (! $file1) {
    die "No old file specified.";
} elsif (! -s $file1) {
    die "$file1 not found or empty.";
} elsif (! $file2) {
    die "No new file specified.";
} elsif (! -s $file2) {
    die "$file2 not found or empty.";
}
# Write the header.
print join("\t", qw(id old-coarse old-fine old-complt old-contam old-good new-coarse new-fine new-complt new-contam new-good)) . "\n";
# This hash will contain the [skcoarse, skfine, complt, contam, good] for each genome in the old file.
my %old;
# Process the old file.
print STDERR "Reading $file1.\n";
open(my $ih, "<$file1") || die "Could not open $file1: $!";
while (! eof $ih) {
    my ($id, $cols) = GetLine($ih);
    $stats->Add(oldLineIn => 1);
    if (defined $id) {
        $old{$id} = $cols;
        $stats->Add(oldLineKept => 1);
    }
}
close $ih; undef $ih;
print STDERR "Reading $file2.\n";
open($ih, "<$file2") || die "Could not open $file2: $!";
while (! eof $ih) {
    my ($id, $cols) = GetLine($ih);
    $stats->Add(newLineIn => 1);
    if (defined $id) {
        $stats->Add(newLineGood => 1);
        my $oldLine = $old{$id};
        if ($oldLine) {
            $stats->Add(newLineMatched => 1);
            print join("\t", $id, @$oldLine, @$cols) . "\n";
        }
    }
}
print STDERR "All done.\n" . $stats->Show();


# Get the data from a file. Return the ID and the list of columns. Return undef for the ID if not all the scores are valid.
sub GetLine {
    my ($ih) = @_;
    my (undef, $id, undef, @fields) = ScriptUtils::get_line($ih);
    my $good = pop @fields;
    my @scores = @fields[8..11];
    if (grep { $_ eq '' } @scores) {
        undef $id;
    }
    return ($id, [@scores, $good]);
}
