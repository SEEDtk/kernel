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
use GPUtils;

=head1 Extract PEG Data From Genome Packages

    good_pegs.pl [ options ] packageDir

This script accepts as input a list of genome IDs and produces an output table of protein-encoding genes from those
genomes. For each PEG, the output file will contain the genome ID, the feature ID, a location string, and a protein
family ID. Optionally, the output can be a FASTA file with the feature ID as the sequence ID and the location string
followed by the protein family ID as the comment.

=head2 Parameters

The only positional parameter is the name of the directory containing the genome packages.

The standard input should be a tab-delimited file.

The command-line options are those found in L<ScriptUtils/ih_options> (to select the standard input file) plus the following.

=over 4

=item col

The index (1-based) of the input colum containing the genome IDs. The default is C<0>, indicating the last column.

=item fasta

If specified, the output file will be in FASTA format rather than tab-delimited format and will include protein sequences.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir', ScriptUtils::ih_options(),
        ['col=i', 'input column index (1-based)', { default => 0 }],
        ['fasta', 'produce FASTA output']
        );
# Get the genome packages directory.
my ($packageDir) = @ARGV;
if (! $packageDir) {
    die "No input packages directory specified.";
} elsif (! -d $packageDir) {
    die "Invalid or missing packages directory $packageDir.";
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the options.
my $col = $opt->col;
my $fastaFormat = $opt->fasta;
# Get all the packages.
my $gHash = GPUtils::get_all($packageDir);
# Loop through the input.
while (! eof $ih) {
    my $genomeID = ScriptUtils::read_col($ih, $col);
    # Get the GTO for the genome.
    my $gto = GPUtils::gto_of($gHash, $genomeID);
    # Extract its pegs.
    my $pegList = GPUtils::all_pegs($gto);
    # Loop through the pegs, producing output.
    for my $peg (@$pegList) {
        # Get the data we need.
        my $fid = $peg->{id};
        my @locs = map { "$_->[0]_$_->[1]$_->[2]$_->[3]" } @{$peg->{location}};
        my $loc = join(",", @locs);
        my $familyData = $peg->{family_assignments};
        if ($familyData) {
            my ($global) = grep { $_->[0] eq 'PGFAM' } @$familyData;
            if ($global) {
                $familyData = $global->[1];
            } else {
                $familyData = '';
            }
        } else {
            $familyData = '';
        }
        # Output according to the format.
        if ($fastaFormat) {
            my $aa = $peg->{protein_translation};
            print ">$fid $loc $familyData\n$aa\n";
        } else {
            print join("\t", $genomeID, $fid, $loc, $familyData) . "\n";
        }
    }
}
