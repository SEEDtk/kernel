package ClosestReps;

use strict;
use warnings;
use Data::Dumper;
use SeedUtils;

sub load_kmers
{
    my ( $roles, $g_to_gs, $k, $shrub ) = @_;

    my $kmer_hash = {};
    my $last_fid = "";
    for my $r (@$roles) {
        # Find all features with this role as their function.
        my @tuples = $shrub->GetAll('Genome Feature Protein AND Feature Feature2Function',
                "Feature2Function(to-link) = ? ORDER BY Feature2Function(from-link), Feature2Function(security) DESC",
                [$r], 'Genome(id) Feature(id) Protein(sequence)');
        for my $tuple (@tuples) {
            my ($g, $fid, $translation) = @$tuple;
            if ($g_to_gs->{$g} && $fid ne $last_fid) {
                my $last = length($translation) - $k;
                for ( my $i = 0 ; ( $i <= $last ) ; $i++ )
                {
                    my $kmer = uc substr( $translation, $i, $k );
                    push( @{ $kmer_hash->{$kmer} }, $g );
                }
                $last_fid = $fid;
            }
        }
    }
    return $kmer_hash;
}

sub process_contigs {
    my($contigs,$kmer_hash,$k) = @_;

    my %hits;
    &add_hits( $contigs, $kmer_hash, $k, \%hits );
    print STDERR "Processing response.\n";
    my $response = &response( \%hits );
    return $response;
}

sub add_hits
{
    my ( $contigs, $kmer_hash, $k, $hits ) = @_;

    my $code = &SeedUtils::standard_genetic_code;
    my $count = scalar @$contigs;
    my $counter = 0;
    foreach my $contig (@$contigs)
    {
        $counter++;
        print STDERR "Processing $contig->[0] + strand ($counter of $count).\n";
        my ( $id, undef, $seq ) = @$contig;
        $seq = uc $seq;
        &add_hits1( $seq, $kmer_hash, $k, $hits, $code );
        print STDERR "Processing $contig->[0] - strand.\n";
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
########################
        my $last1 = ($k-1) * 3;
        my $prot = 'x' x $k;
        my $j;
        for ($j=0; ($j <= $last1); $j += 3)
        {
            my $dna = uc substr($seq,$i+$j,3);
            my $aa = $code->{$dna};
            if (! $aa) { $aa = 'x' }
            substr($prot,$j/3,1) = $aa;
        }
#
# I tried to speed things up by replacing the following two lines
# with the above.  Not much help.
#
#        my $dna = substr( $seq, $i, $k * 3 );
#        my $prot = &SeedUtils::translate( $dna, $code, 0 );
########################

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
