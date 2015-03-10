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
use Data::Dumper;
use SeedUtils;

=head1 project reference genome to a close strain

    fast_project -r RererenceGenomeDir -g SkeletalGenomeDir [ options ]

project a reference genome to call features


=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -r ReferenceGenomeDir

a path to a SEED genome directory for the reference genome

=item -g SkeletalGenomeDir

a path to a skeletal SEED genome directory that must include

    contigs
    GENETIC_CODE (if not 11)

=back

=cut

my $k = 30;    # kmers for generating map.  Chosen conservatively.

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'reference|r=s', 'Path to Reference Genome Directory' ],
    [ 'genome|g=s',    'Path to a skeletal genome directory' ]
);
my $refD    = $opt->reference;
my $genomeD = $opt->genome;

my @ref_tuples = &gjo::seqlib::read_fasta("$refD/contigs");
my @g_tuples   = &gjo::seqlib::read_fasta("$genomeD/contigs");

my $map = &build_mapping( \@ref_tuples, \@g_tuples );
&nuild_features( $map, $refD, $genomeD );

sub build_mapping
{
    my ( $r_contigs, $g_contigs ) = @_;

    my $g_hash = &build_g_hash( $g_contigs, $k );
    my $pins = &build_pins( $r_contigs, $k, $g_hash );
    my $map = &fill_pins( $pins, \@ref_tuples, \@g_tuples );

    return $map;
}

sub build_g_hash
{
    my ( $g_contigs, $k );

    my $g_hash = {};
    my %seen;
    foreach my $tuple (@$g_contigs)
    {
        my ( $contig_id, $comment, $seq ) = @$tuple;
        my $last = length($seq) - $k;
        for ( my $i = 0 ; ( $i <= $last ) ; $i++ )
        {
            my $kmer = uc substr( $seq, $i, $k );
            if ( $kmer !~ /^[ACGT]/ )
            {
                my $comp = &rev_comp($kmer);
                if ( $g_hash->{$kmer} )
                {
                    $seen{$kmer} = 1;
                    $seen{$comp} = 1;
                }
                $g_hash->{$kmer} = [ $contig_id, $i, "+" ];
                $g_hash->{$comp} = [ $contig_id, $i + $k - 1, "-" ];
            }
        }
    }

    foreach my $kmer ( keys(%seen) )
    {
        delete $g_hash->{$kmer};
        delete $g_hash->{$kmer};
    }
    return $g_hash;
}

sub build_pins
{
    my ( $r_contigs, $k, $g_hash ) = @_;

    my @pins;
    foreach my $tuple (@$r_contigs)
    {
        my ( $contig_id, $comment, $seq ) = @$tuple;
        my $last = length($seq) - $k;
        for ( my $i = 0 ; ( $i <= $last ) ; $i++ )
        {
            my $kmer = uc substr( $seq, $i, $k );
            if ( $kmer !~ /^[ACGT]/ )
            {
                my $g_pos = $g_hash->{$kmer};
                if ($g_pos)
                {
                    push( @pins, [ [ $contig_id, $i, '+' ], $g_pos ] );
                }
            }
        }
    }
    @pins = sort { ( $a->[0] cmp $b->[0] ) or ( $a->[1] <=> $b->[1] ) } @pins;
    return \@pins;
}

sub fill_pins
{
    my ( $pins, $ref_tuples, $g_tuples ) = @_;

    my %ref_seqs = map { ( $_->[0] => $_->[2] ) } @$ref_tuples;
    my %g_seqs   = map { ( $_->[0] => $_->[2] ) } @$g_tuples;

    my @filled;
    for ( my $i = 0 ; ( $i < @$pins ) ; $i++ )
    {
        if ( $i == @{$pins} )
        {
            push( @filled, $pins->[$i] );
        }
        else
        {
            my @expanded = &fill_between( $pins->[$i], $pins->[ $i + 1 ],
                \%ref_seqs, \%g_seqs );
            push( @filled, @expanded );
        }
    }
    return @filled;
}

sub fill_between
{
    my ( $pin1, $pin2, $ref_seqs, $g_seqs ) = @_;

    my ( $rp1, $gp1 ) = @$pin1;
    my ( $rp2, $gp2 ) = @$pin2;
    my ( $contig_r_1, $pos_r_1, $strand_r_1 ) = @$rp1;
    my ( $contig_r_2, $pos_r_2, $strand_r_2 ) = @$rp2;
    my ( $contig_g_1, $pos_g_1, $strand_g_1 ) = @$gp1;
    my ( $contig_g_2, $pos_g_2, $strand_g_2 ) = @$gp2;

    my @expanded = ($pin1);
    if (
           ( $contig_r_1 eq $contig_r_2 )
        && ( $contig_g_1 eq $contig_g_2 )
        && ( $strand_g_1 eq $strand_g_2 )
        && ( ( $pos_r_2 - $pos_r_1 ) == abs( $pos_g_2 - $pos_g_1 ) )
        && ( ( $pos_r_2 - $pos_r_1 ) > 1 )
        && &same(
            [ $contig_r_1, '+', $pos_r_1, $pos_r_2 - 1, $ref_seqs ],
            [
                $contig_g_1,
                $strand_g_1,
                ( $strand_g_1 eq '+' )
                ? ( $pos_g_1, $pos_g_2 - 1 )
                : ( $pos_g_1, $pos_g_2 + 1 )
            ]
        )
      )
    {
        my $p_r = $pos_r_1;
        my $p_g = $pos_g_1;
        while ( $p_r < $pos_r_2 )
        {
            push(
                @expanded,
                [
                    [ $contig_r_1, '+',         $p_r ],
                    [ $contig_g_1, $strand_g_1, $p_g ]
                ]
            );
            $p_r++;
            $p_g = ( $strand_g_1 eq "+" ) ? $p_g + 1 : $p_g - 1;
        }
    }
    return @expanded;
}

sub same
{
}
