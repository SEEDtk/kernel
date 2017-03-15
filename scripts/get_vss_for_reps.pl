use strict;
use Data::Dumper;
use gjoseqlib;
use SeedUtils;
use RepKmers;

my $usage = "usage: get_vss_for_reps Data\n";
my $dataD;
(
 ($dataD = shift @ARGV)
)
    || die $usage;

my $K;
if (open(K,"<$dataD/K"))
{
    $K = <K>;
    chomp $K
}
else
{
    $K = 8;
}

my %tails;
open(FASTA,"<$dataD/6.1.1.20.fasta") 
    || die "could not open $dataD/6.1.1.20.fasta";
my $entry;
my %seqs;

while ($entry = &gjoseqlib::read_next_fasta_seq(\*FASTA))
{
    my($id,undef,$seq) = @$entry;
    $seqs{$id} = $seq;
    my $tail = lc substr($seq,-12);
    push(@{$tails{$tail}},$id);
}
close(FASTA);
print STDERR "loaded tails\n";
open(VSS,">$dataD/vss") || die "could not open $dataD/vss";
my %seen;

foreach my $tail (keys(%tails))
{
    my $hits = $tails{$tail};
    if (@$hits > 1)
    {
        $_ = @$hits; print STDERR "sz hit = $_\n";
	my $longest = &longest($hits,\%seqs);
	if (! $seen{$longest})
	{
	    my $seq1 = $seqs{$longest};
	    my @vss = $longest;
	    foreach my $id (keys(%seqs))
	    {
		if (($id ne $longest) && (! $seen{$id}))
		{
		    my $common = &RepKmers::sim($seq1,$seqs{$id},$K);
		    if ($common >= 200)
		    {
			push(@vss,$id);
		    }
		}
	    }
	    if (@vss > 1)
	    {
		$_ = @vss; print STDERR "vss=$_\n";
		foreach $_ (sort @vss)
		{
		    print VSS $_,"\n";
		    $seen{$_} = 1;
		}
		print VSS "//\n";
	    }
	}
    }
}
close(VSS);

sub longest {
    my($hits,$seqs) = @_;

    my $sofar;
    my $best;
    foreach my $hit (@$hits)
    {
	my $len = length ($seqs->{$hit});
	if ((! $best) || ($len >= $sofar))
	{
	    $best = $hit;
	    $sofar = $len;
	}
    }
    return $best;
}
