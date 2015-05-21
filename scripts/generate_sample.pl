use strict;
use Data::Dumper;
use List::Util qw(shuffle);
use Shrub::Contigs;
use Shrub;
my $shrub = new Shrub;

my $Ng;
my $Na;

(
 ($Ng = shift @ARGV) &&
 ($Na = shift @ARGV)
)
    || die "usage: generate_sample NumberGenomes MaxAbundance";
my @genomes = map  { ($_ =~ /^(\S+)/) ? $1 : () } &shuffle(`all_genomes`);
$#genomes = $Ng-1;

open(LOG,">log") || die "bad";
my $id = 1;
foreach my $g (@genomes)
{
    my $contigO = Shrub::Contigs->new($shrub,$g);
    my @triples = $contigO->tuples;
    my $real_ln = 0;
    foreach $_ (@triples)
    {
	$real_ln += length($_->[2]);
    }
    my $cov = int(rand() * $Na);
    print LOG "genome=$g\ncov=$cov\n";
    my $out_ln = 0;
    while ($out_ln < ($cov * $real_ln))
    {
	my $desired = int(rand() * 10000) + 200;
	my $start = int(rand() * scalar @triples);
	my $got = 0;
	while (! $got)
	{
	    my $len_of_contig = length($triples[$start]->[2]);
	    if ($desired > $len_of_contig)
	    {
		$start++;
		if ($start == @triples) { $start = 0 }
	    }
	    else
	    {
		my $first = int(rand() * ($len_of_contig - $desired));
		my $seq = substr($triples[$start]->[2],$first,$desired);
		print ">",$id++,"_length_",$desired,"_covg_",$cov,"\n";
		print $seq,"\n";
		$out_ln += $desired;
		$got = 1;
	    }
	}
    }
}
