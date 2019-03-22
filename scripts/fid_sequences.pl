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
use Shrub;
use ScriptUtils;
use Shrub::Contigs;

=head1 Add DNA and Protein Sequences to a Feature File

    fid_sequences.pl [ options ]

This script takes as input a tab-delimited file containing feature IDs and appends DNA and/or protein sequences.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> (to specify the L<Shrub> database) and
L<ScriptUtils/ih_options> (to specify the standard input) plus the following.

=over 4

=item col

The index (1-based) of the input column containing the feature IDs.  The default is C<0>, indicating the last column.

=item prot

If specified, protein sequences will be output.  This is the default if neither C<--prot> or C<--dna>.

=item dna

If specified, DNA sequences will be output.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(), ScriptUtils::ih_options(),
        ['col|c=i', 'index (1-based) of column with feature IDs', { default => 0 }],
        ['prot|p', 'output protein sequences'],
        ['dna|n', 'output DNA sequences']
        );
# Compute the outputs.
my $dna = $opt->dna;
my $prot = $opt->prot;
# Default to protein only.
if (! $dna && ! $prot) {
    $prot = 1;
}
# Compute the feature ID column.
my $fidCol = $opt->col - 1;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# This hash will map each genome ID to a list of couplets.  A couplet in this case consists of [fid,line].
my %genomes;
# Read in the whole input and sort it by genome ID.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    my @flds = split /\t/, $line;
    my $fid = $flds[$fidCol];
    if ($fid =~ /^fig\|(\d+\.\d+)/) {
        push @{$genomes{$1}}, [$fid, $line];
    }
}
# Now loop through the genomes, processing each one individually.
for my $genome (sort keys %genomes) {
    # Get this genome's couplets.
    my $couplets = $genomes{$genome};
    # Extract the features.
    my %fids = map { $_->[0] => [] } @$couplets;
    if ($dna) {
        # Get the DNA for these features.
        my $contigs = Shrub::Contigs->new($shrub, $genome);
        for my $fid (keys %fids) {
            $fids{$fid} = [$contigs->fdna($fid)];
        }
    }
    if ($prot) {
        # Get the proteins for these features.
        for my $fid (keys %fids) {
            my ($pSeq) = $shrub->GetFlat('Feature2Protein Protein', 'Feature2Protein(from-link) = ?', [$fid], 'Protein(sequence)');
            $pSeq //= '';
            push @{$fids{$fid}}, $pSeq;
        }
    }
    # Process the couplets again, writing output.
    for my $couplet (@$couplets) {
        my ($fid, $line) = @$couplet;
        my $seqsL = $fids{$fid};
        print join("\t", $line, @$seqsL) . "\n";
    }
}

## TODO process the input to produce the output.