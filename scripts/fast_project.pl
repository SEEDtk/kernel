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
use Shrub;
use ScriptUtils;
use gjo::seqlib;

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

my $k = 4;    # kmers for generating map.  Chosen conservatively.

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

my $genetic_code = &get_genetic_code($genomeD);

my @ref_tuples = &gjo::seqlib::read_fasta("$refD/contigs");
my @g_tuples   = &gjo::seqlib::read_fasta("$genomeD/contigs");

my $map = &build_mapping( \@ref_tuples, \@g_tuples );
&build_features( $map, $refD, $genomeD, \@g_tuples, $genetic_code );

sub build_mapping {
	my ( $r_contigs, $g_contigs ) = @_;
	print STDERR "Calling build_hash for reference\n";
	my $r_hash = &build_hash( $r_contigs, $k );
	print STDERR "Calling build_hash for new\n";
	my $g_hash = &build_hash( $g_contigs, $k );

	my $pins = &build_pins( $r_contigs, $k, $g_hash, $r_hash );
	my @map = &fill_pins( $pins, \@ref_tuples, \@g_tuples );

	return \@map;
}

# a hash has a 0-base for each kmer (kmer is a key to a 0-based location)
sub build_hash {
	my ( $contigs, $k ) = @_;

	my $hash = {};
	my %seen;
	foreach my $tuple (@$contigs) {
		my ( $contig_id, $comment, $seq ) = @$tuple;
		my $last = length($seq) - $k;
		for ( my $i = 0 ; ( $i <= $last ) ; $i++ ) {
			my $kmer = uc substr( $seq, $i, $k );
			if ( $kmer !~ /[^ACGT]/ ) {
				my $comp = &rev_comp($kmer);
				if ( $hash->{$kmer} ) {
					$seen{$kmer} = 1;
					$seen{$comp} = 1;
				}
				$hash->{$kmer} = [ $contig_id, "+", $i ];
				$hash->{$comp} = [ $contig_id, "-", $i + $k - 1 ];
			}
		}
	}

	foreach my $kmer ( keys(%seen) ) {
		delete $hash->{$kmer};
	}
	print STDERR &Dumper('hash',$hash);
	return $hash;
}

# pins are 0-based 2-tuples
sub build_pins {
	my ( $r_contigs, $k, $g_hash, $r_hash ) = @_;

	my @pins;
	foreach my $tuple (@$r_contigs) {
		my ( $contig_id, $comment, $seq ) = @$tuple;
		my $last = length($seq) - $k;

		my $i = 0;
		while ( $i <= $last ) {
			my $kmer = uc substr( $seq, $i, $k );
			if ( ( $kmer !~ /[^ACGT]/ ) && $r_hash->{$kmer} ) {
				my $g_pos = $g_hash->{$kmer};
				if ($g_pos) {
					my ( $g_contig, $g_strand, $g_off ) = @$g_pos;
					for ( my $j = 0 ; $j < $k ; $j++ ) {
						if ( $g_strand eq '+' ) {
							push(
								@pins,
								[
									[ $contig_id, '+', $i + $j ],
									[ $g_contig,  '+', $g_off + $j ]
								]
							);
						}
						else {
							push(
								@pins,
								[
									[ $contig_id, '+', $i + $j ],
									[ $g_contig,  '-', $g_off - $j ]
								]
							);
						}
					}
					$i = $i + $k;
				}
				else {
					$i++;
				}
			} else { $i++}
		}
	}
	print STDERR &Dumper(['0-based pins', \@pins] );
	@pins = sort {
		     ( $a->[0]->[0] cmp $b->[0]->[0] )
		  or ( $a->[0]->[2] <=> $b->[0]->[2] )
	} @pins;
	return \@pins;
}

sub fill_pins {
	my ( $pins, $ref_tuples, $g_tuples ) = @_;

	my %ref_seqs = map { ( $_->[0] => $_->[2] ) } @$ref_tuples;
	my %g_seqs   = map { ( $_->[0] => $_->[2] ) } @$g_tuples;

	my @filled;
	for ( my $i = 0 ; ( $i < @$pins ) ; $i++ ) {
		if ( $i == ( @$pins - 1 ) ) {
			push( @filled, $pins->[$i] );
		}
		else {
			my @expanded = &fill_between( $pins->[$i], $pins->[ $i + 1 ],
				\%ref_seqs, \%g_seqs );
			push( @filled, @expanded );
		}
	}
	return @filled;
}

sub fill_between {
	my ( $pin1, $pin2, $ref_seqs, $g_seqs ) = @_;

	my ( $rp1, $gp1 ) = @$pin1;
	my ( $rp2, $gp2 ) = @$pin2;
	my ( $contig_r_1, $strand_r_1, $pos_r_1 ) = @$rp1;
	my ( $contig_r_2, $strand_r_2, $pos_r_2 ) = @$rp2;
	my ( $contig_g_1, $strand_g_1, $pos_g_1 ) = @$gp1;
	my ( $contig_g_2, $strand_g_2, $pos_g_2 ) = @$gp2;

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
				: ( $pos_g_1, $pos_g_2 + 1 ),
				$g_seqs
			]
		)
	  )
	{
		my $p_r = $pos_r_1;
		my $p_g = $pos_g_1;
		while ( $p_r < $pos_r_2 ) {
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

sub same {
	my ( $gap1, $gap2 ) = @_;
	my ( $c1, $b1, $e1, $seqs1 ) = @$gap1;
	my ( $c2, $b2, $e2, $seqs2 ) = @$gap2;

	my $seq1 = &seq_of( $c1, $b1, $e1, $seqs1 );
	my $seq2 = &seq_of( $c2, $b2, $e2, $seqs2 );
	if ( length($seq1) < 20 ) {
		return 1;
	}
	else {
		my $iden = 0;
		my $len  = length($seq1);
		for ( my $i = 0 ; ( $i < $len ) ; $i++ ) {
			if ( substr( $seq1, $i, 1 ) eq substr( $seq2, $i, 1 ) ) {
				$iden++;
			}
		}
		return ( ( $iden / $len ) >= 0.8 );
	}
}

sub seq_of {
	my ( $c, $b, $e, $seqs ) = @_;

	my $seq = $seqs->{$c};
	if ( $b <= $e ) {
		return uc substr( $seq, $b, ( $e - $b ) + 1 );
	}
	else {
		return uc &rev_comp( substr( $seq, $e, ( $b - $e ) + 1 ) );
	}
}

sub build_features {
	my ( $map, $refD, $genomeD, $g_tuples, $genetic_code ) = @_;

	my %g_seqs = map { ( $_->[0] => $_->[2] ) } @$g_tuples;

	my %refH;
	foreach my $pin (@$map) {
		my ( $ref_loc, $g_loc ) = @$pin;
		my ( $r_contig, $r_strand, $r_pos ) = @$ref_loc;
		$refH{ $r_contig . ",$r_pos" } = $g_loc;
	}
	print STDERR "refH in build features", Dumper( \%refH );
	mkdir( "$genomeD/Features", 0777 );
	opendir( FEATURES, "$refD/Features" )
	  || die "could not open $refD/Features";
	my @types = grep { $_ !~ /^\./ } readdir(FEATURES);
	closedir(FEATURES);

	foreach my $type (@types) {
		my $dir = "$genomeD/Features/$type";

		my %deleted_features;
		if ( -s "$refD/Features/$type/deleted.features" ) {
			%deleted_features =
			  map { ( $_ =~ /^(\S+)/ ) ? ( $1 => 1 ) : () }
			  `cat $refD/Features/$type/deleted.features`;
		}
		mkdir( $dir, 0777 );
		if (   open( TBL, ">$dir/tbl" )
			&& open( FASTA, ">$dir/fasta" ) )
		{
			foreach my $f_line (`cat $refD/Features/$type/tbl`) {
				print STDERR $f_line;
				if (   ( $f_line =~ /^(\S+)\t(\S+)_(\S+)_(\S+)\t/ )
					&& ( !$deleted_features{$1} ) )
				{
					my ( $fid, $r_contig, $r_beg, $r_end ) =
					  ( $1, $2, $3 - 1, $4 - 1 );
					if (   ( my $g_locB = $refH{ $r_contig . ",$r_beg" } )
						&& ( my $g_locE = $refH{ $r_contig . ",$r_end" } ) )
					{
						my ( $g_contig1, $g_strand1, $g_pos1 ) = @$g_locB;
						my ( $g_contig2, $g_strand2, $g_pos2 ) = @$g_locE;
						if (
							   ( $g_contig1 eq $g_contig2 )
							&& ( $g_strand1 eq $g_strand2 )
							&& (
								abs( $g_pos1 - $g_pos2 ) ==
								abs( $r_beg - $r_end ) )
						  )
						{
							my $g_location = join( "_",
								( $g_contig1, $g_pos1 + 1, $g_pos2 + 1 ) );
							print STDERR "Adding $g_location\n";
							my $seq =
							  &seq_of_feature( $type, $genetic_code, $g_contig1,
								$g_pos1, $g_pos2, \%g_seqs );
							print STDERR $seq, "\n";
							if ($seq) {
								print TBL "$fid\t$g_location\t\n";
								$r_beg++;
								$r_end++;
								print FASTA
								  ">$fid $r_contig $r_beg $r_end\n$seq\n";
							}
						}
					}
				}
			}
		}
		close(TBL);
		close(FASTA);
	}
}

sub get_genetic_code {
	my ($dir) = @_;

	if ( !-s "$dir/GENETIC_CODE" ) { return 11 }
	my @tmp = `cat $dir/GENETIC_CODE`;
	chomp $tmp[0];
	return $tmp[0];
}

sub seq_of_feature {
	my ( $type, $genetic_code, $g_contig, $g_beg, $g_end, $g_seqs ) = @_;
	my $dna = &seq_of( $g_contig, $g_beg, $g_end, $g_seqs );
	if ( ( $type ne "peg" ) && ( $type ne "CDS" ) ) {
		return $dna;
	}
	else {
		my $code = &SeedUtils::standard_genetic_code;
		if ( $genetic_code == 4 ) {
			$code->{"TGA"} = "W";    # code 4 has TGA encoding tryptophan
		}
		my $tran = &SeedUtils::translate( $dna, $code, 1 );
		$tran =~ s/\*$//;
		return ( $tran =~ /\*/ ) ? undef : $tran;
	}
}
