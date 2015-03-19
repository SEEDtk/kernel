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
use Shrub;
use ScriptUtils;
use Data::Dumper;
use gjo::seqlib;

=head1 Get closest genomes in coreSEED

    get_closest_coreSEED_genomes [ options ] < requests

guts of a server that processes requests for closest coreSEED genomes
to new genomes.

The input will be requests of the form:
    
    contigs-in-fasta 
    //
    next-set-of-contigs-in-fasta
    //
    .
    .
    .

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=cut

my $k = 10;

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '', Shrub::script_options(), ScriptUtils::ih_options(),
    [] );
my $ih = ScriptUtils::IH( $opt->input );

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my @functions = <DATA>;
chomp @functions;
my $kmer_hash = &load_kmers( \@functions, $shrub, $k );

while ( my $contigs = &read_contigs )
{
    my %hits;
    &add_hits( $contigs, $kmer_hash, $k, \%hits );
    my $response = &response( \%hits );
    print $response, "\n//\n";
}

sub load_kmers
{
    my ( $functions, $shrub, $k ) = @_;
    $kmer_hash = {};
    foreach my $function (@$functions)
    {
        my @tuples = $shrub->GetAll(
            "Function Function2Feature Feature Feature2Protein Protein",
            "(Function(description) = ?) AND (Function2Feature(security) = ?)",
            [ $function, 2 ],
            "Feature(id) Protein(sequence)"
        );
        foreach my $tuple (@tuples)
        {
            my ( $peg, $translation ) = @$tuple;
            my $g    = &SeedUtils::genome_of($peg);
            my $last = length($translation) - $k;
            for ( my $i = 0 ; ( $i <= $last ) ; $i++ )
            {
                my $kmer = uc substr( $translation, $i, $k );
                push( @{ $kmer_hash->{$kmer} }, $g );
            }
        }
    }
    return $kmer_hash;
}

sub read_contigs
{
    my $contigs = [];
    $/ = "\n//\n";
    my $req = <$ih>;
    if ($req)
    {
        chomp $req;
        my @entries = split( /\n\>\s*/, $req );
        foreach my $entry (@entries)
        {
            if ( $entry =~ /^(>\s*)?(\S+)[^\n]*\n(.*)/s )
            {
                my ( $id, $seq ) = ( $2, $3 );
                $seq =~ s/\s+//gs;
                push( @$contigs, [ $id, '', $seq ] );
            }
        }
    }
    $/ = "\n";
    if (@$contigs)
    {
        return $contigs;
    }
    return undef;
}

sub add_hits
{
    my ( $contigs, $kmer_hash, $k, $hits ) = @_;

    my $code = &SeedUtils::standard_genetic_code;
    foreach my $contig (@$contigs)
    {
        my ( $id, undef, $seq ) = @$contig;
        $seq = uc $seq;
        &add_hits1( $seq, $kmer_hash, $k, $hits, $code );
        my $seqR = &SeedUtils::rev_comp($seq);
        &add_hits1( $seqR, $kmer_hash, $k, $hits, $code );
    }
}

sub add_hits1
{
    my ( $seq, $kmer_hash, $k, $hits, $code ) = @_;

    my $i = 0;
    my $last = length($seq) - ( 3 * $k );
    for ( $i = 0 ; ( $i <= $last ) ; $i++ )
    {
        my $dna = substr( $seq, $i, $k * 3 );
        my $prot = &SeedUtils::translate( $dna, $code, 0 );
        my $occ = $kmer_hash->{$prot};
        if ($occ)
        {
            foreach $_ (@$occ)
            {
                $hits->{$_}++;
            }
        }
    }
}

sub response
{
    my ($hits) = @_;

    my @poss =
      grep { $hits->{$_} >= 5 }
      sort { $hits->{$b} <=> $hits->{$a} } keys(%$hits);
    if ( @poss == 0 ) { return 'none' }
    return join( "\n", map { "$hits->{$_}, $_" } @poss);
}

__DATA__
Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
Phenylalanyl-tRNA synthetase beta chain (EC 6.1.1.20)
Preprotein translocase secY subunit (TC 3.A.5.1.1)
GTP-binding and nucleic acid-binding protein YchF
Translation initiation factor 2
Signal recognition particle, subunit Ffh SRP54 (TC 3.A.5.1.1)
Histidyl-tRNA synthetase (EC 6.1.1.21)
Methionyl-tRNA synthetase (EC 6.1.1.10)
Isoleucyl-tRNA synthetase (EC 6.1.1.5)
Valyl-tRNA synthetase (EC 6.1.1.9)
Seryl-tRNA synthetase (EC 6.1.1.11)
Alanyl-tRNA synthetase (EC 6.1.1.7)
