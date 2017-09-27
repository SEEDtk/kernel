use strict;
use Data::Dumper;

# usage: collapse < shared
#
#                   where shared is made up of lines
#                   containing
#
#        [genome1,genome2,kmers-in-common]

my %shared;
my $ih;
if ($ARGV[0]) {
    open($ih, "<$ARGV[0]") || die "Could not open input: $!";
} else {
    $ih = \*STDIN;
}
foreach $_ (<$ih>)
{
    chop;
    my($g1,$g2,$in_common) = split(/\t/,$_);
    $shared{$g1}->{$g2} = $shared{$g2}->{$g1} = $in_common;
}

my @nodes = map { [$_] } sort keys(%shared);

while (@nodes > 1)
{
    my $best_pair;
    my $bestI;
    my $bestJ;
    my($i,$j);
    for ($i=0; ($i < (@nodes - 1)); $i++)
    {
        for ($j=$i+1; ($j< @nodes); $j++)
        {
            my $avg_in_common = &avg_in_common(\%shared,$nodes[$i],$nodes[$j]);
            if ((! $best_pair) ||
                ($avg_in_common->{COMMON} > $best_pair->{COMMON}))
            {
                $best_pair = $avg_in_common;
                $bestI = $i;
                $bestJ = $j;
            }
        }
    }
    my $genomes1 = $nodes[$bestI];
    my $genomes2 = $nodes[$bestJ];
    my @new = sort (@$genomes1,@$genomes2);
    &display($genomes1,$genomes2,$best_pair);
    splice(@nodes,$bestJ,1);
    $nodes[$bestI] = [@new];
}

sub display {
    my($genomes1,$genomes2,$best_pair) = @_;

    print "merging\n",join(",",@$genomes1), "\n",
                      join(",",@$genomes2), "\n";
    print "average in common: ",$best_pair->{AVG},"\n\n";
}

sub avg_in_common {
    my($shared,$n1,$n2) = @_;

    my $n = 0;
    my $tot = 0;
    foreach my $x (@$n1)
    {
        if ($shared->{$x})
        {
            foreach my $y (@$n2)
            {
                if ($shared->{$x}->{$y})
                {
                    $n++;
                    $tot += $shared->{$x}->{$y};
                }
            }
        }
    }
    return { COMMON => $tot,
             AVG    => sprintf("%0.2f",$tot/$n)
    }
}

