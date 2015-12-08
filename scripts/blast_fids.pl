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
use gjoseqlib;
use gjo::BlastInterface;
use SeedUtils;
use Shrub::Contigs;
use BasicLocation;

=head1 BLAST Features Against Genomes or Themselves

    blast_fids.pl [ options ] genomeID

This script blasts features against other features. The query features are specified as a list of feature IDs.
The database features can be the same as the query features or they can be taken from a single genome. Either
C<blastn> or C<blastp> can be selected, meaning that we are always BLASTing proteins against proteins or DNA
against DNA.

=head2 Parameters

The single positional parameter is a genome ID or FASTA file name. If a genome ID is specified, then the genome is
used as the database for the BLAST. If a FASTA file name is specified, the FASTA file is used. Otherwise, the input
sequences are blasted against themselves.

The command-line parameters are those specified in L<Shrub/script_options> (database connection) and L<ScriptUtils/ih_options>
(standard input) plus the following.

=over 4

=item dna

If specified, then DNA blasting is used. This parameter is mutually exclusive with C<prot>.

=item prot

If specified, then protein blasting is used. This parameter is the default, and is mutually exclusive with C<dna>.

=item fasta

If specified, then the list of features is specified as a FASTA file. This parameter is mutually exclusive with C<col>.
The FASTA file will be on the standard input.

=item col

The column from the standard input from which the query sequence feature IDs is to be taken. The default is the last
column.

=item hsp

If specified, then the output is in the form of HSP data (see L<Hsp>). This is the default, and is mutually exclusive with C<sim>.

=item sim

If specified, then the output is in the form of similarity data (see L<Sim>). This parameter is mutually exclusive with C<hsp>.

=item BLAST Parameters

The following may be specified as BLAST parameters

=over 8

=item maxE

Maximum E-value (default C<1e-10>).

=item maxHSP

Maximum number of returned results (before filtering). The default is to return all results.

=item minScr

Minimum required bit-score. The default is no minimum.

=item percIdentity

Minimum percent identity. The default is no minimum.

=item minLen

Minimum permissible match length (used to filter the results). The default is no filtering.

=back

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('genomeID',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['mode' => hidden => { one_of => [ ['dna|n' => 'use DNA blasting'], ['prot|p' => 'use protein blasting'] ] }],
        ['fasta|f', 'if specified, input is assumed to be a fasta file'],
        ['col|c=i', 'if specified, index (1-based) of input column containing feature IDs'],
        ['output' => hidden => { one_of => [ [ 'hsp' => 'produce HSP output'], ['sim' => 'produce similarity output'] ]}],
        ['maxE|e', 'maximum e-value', { default => 1e-10 }],
        ['maxHSP|b', 'if specified, the maximum number of returned results (before filtering)'],
        ['minScr=f', 'if specified, the minimum permissible bit score'],
        ['percIdentity=f', 'if specified, the minimum permissible percent identity'],
        ['minLen|l=i', 'if specified, the minimum permissible match lengt (for filtering)']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Compute the blast mode.
my $mode = $opt->mode // 'prot';
# This variable will contain the blast program to use.
my $blastProg = ($mode eq 'dna' ? 'blastn' : 'blastp');
# This hash contains the BLAST parameters.
my %blast;
$blast{outForm} = $opt->output // 'hsp';
$blast{maxE} = $opt->maxe;
$blast{maxHSP} = $opt->maxhsp // 0;
$blast{minIden} = $opt->percidentity // 0;
$blast{minLen} = $opt->minlen // 0;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# We need the query features. These can be specified as a list of feature IDs from the standard input, or as
# a FASTA file. The query sequences will be put into the following variable as a list of triples.
my @query;
if ($opt->fasta) {
    # Here the query features are input as a FASTA file.
    @query = gjoseqlib::read_fasta($ih);
} else {
    # Here we are expecting feature IDs.
    my $col = $opt->col;
    my @fids;
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        my @flds = split /\t/, $line;
        if ($col) {
            push @fids, $flds[$col - 1];
        } else {
            push @fids, pop @flds;
        }
    }
    # What we do next depends on whether we are processing DNA or proteins.
    if ($mode eq 'prot') {
        # Proteins are easy. Just get the sequence for each feature.
        for my $fid (@fids) {
            my ($seq) = $shrub->GetFlat('Feature2Protein Protein', 'Feature2Protein(from-link) = ?',
                    [$fid], 'Protein(sequence)');
            if ($seq) {
                push @query, [$fid, '', $seq];
            } else {
                print STDERR "No protein sequence found for $fid.\n";
            }
        }
    } else {
        # For DNA, we need to extract the DNA from the contigs. For performance, we sort the feature IDs
        # by genome.
        my %gFids;
        for my $fid (@fids) {
            my $genome = SeedUtils::genome_of($fid);
            if (! $genome) {
                print STDERR "$fid has an invalid format.\n";
            } else {
                push @{$gFids{$genome}}, $fid;
            }
        }
        # Loop through the genomes found.
        for my $genome (keys %gFids) {
            # Load this genome's DNA.
            my $contigs = Shrub::Contigs->new($shrub, $genome);
            # Loop through the features.
            my $fidList = $gFids{$genome};
            for my $fid (@$fidList) {
                my $seq = $contigs->fdna($fid);
                push @query, [$fid, '', $seq];
            }
        }
    }
}
# Now @query contains an array of FASTA triples for the feature sequences.
# We are done with the standard input.
close $ih;
# The next goal is to get the database.
my $db;
# Check for a positional parameter.
my ($blastdb) = @ARGV;
if (! $blastdb) {
    # No positional parameter, so we are blasting the fids against themselves.
    $db = \@query;
} elsif ($blastdb =~ /^\d+\.\d+$/) {
    # Here the database is a genome. Read in the triples for it.
    my @triples;
    if ($mode eq 'prot') {
        $shrub->write_prot_fasta($blastdb, \@triples);
    } else {
        $shrub->write_peg_fasta($blastdb, \@triples);
    }
    $db = \@triples;
}
# Now run the BLAST.
my $matches = gjo::BlastInterface::blast(\@query, $db, $blastProg, \%blast);
# Format the output.
for my $match (@$matches) {
    print join("\t", @$match) . "\n";
}
