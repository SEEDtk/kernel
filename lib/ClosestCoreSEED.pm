package ClosestCoreSEED;

use strict;
use warnings;
use Data::Dumper;
use SeedUtils;

sub load_kmers
{
    my ( $functions, $shrub, $k ) = @_;
    my $kmer_hash = {};
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

sub process_contigs {
    my($contigs,$kmer_hash,$k) = @_;

    my %hits;
    &add_hits( $contigs, $kmer_hash, $k, \%hits );
    my $response = &response( \%hits );
    return $response;
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

1;
