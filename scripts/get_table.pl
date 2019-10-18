use strict;
use Data::Dumper;
use BlastUtils;
use gjoseqlib;

my %seen;
my @all_seqs    = &gjoseqlib::read_fasta('EvalBinning/6.1.1.20.fasta');
my @seqs;
foreach my $seq (@all_seqs)
{
    if (($seq->[1] =~ /^\d+\.\d+\t(\S+ \S+)/) && (! $seen{$1}))
    {
        $seen{$1} = 1;
        $seq->[1] = $1;
        push(@seqs,$seq);
    }
}
gjoseqlib::write_fasta('EvalBinning/condensed.fa', \@seqs);
#my %genomes = map { ($_->[1] =~ /(\d+\.\d+)\t(\S+ \S+)/) ? ($1 => $2) : () } @seqs;
#foreach my $seq (@seqs)
#{
#    my @matches = grep { ($_->[0] ne $_->[1]) } &BlastUtils::blastp([$seq],\@seqs,{});
#    my $g1;
#    my $g2;
#    my $iden = $matches[0]->[2];
#    if (($matches[0]->[0] =~ /^fig\|(\d+\.\d+)/) && ($g1 = $1) &&
#	($matches[0]->[1] =~ /^fig\|(\d+\.\d+)/) && ($g2 = $1))
#    {
#        my $gn1 = $genomes{$g1};
#	my $gn2 = $genomes{$g2};
#	if ($gn2 lt $gn1) { ($gn1,$gn2) = ($gn2,$gn1) }
#	if ($gn1 ne $gn2)
#	{
#	    my $pair = "$gn1\t$gn2";
#	    if (! $seen{$pair})
#	    {
#		$seen{$pair} = 1;
#		print join("\t",($iden,$g1,$gn1,$g2,$gn2)),"\n";
#	    }
#	}
#    }
#}

