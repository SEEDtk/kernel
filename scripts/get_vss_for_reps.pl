use strict;
use Data::Dumper;
use gjoseqlib;
use SeedUtils;
use RepKmers;
use Stats;

my $stats = Stats->new();
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
print STDERR "Kmer size is $K.\n";

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
    $stats->Add(tailStored => 1);
}
close(FASTA);
my $tailCount = scalar keys %tails;
print STDERR "loaded tails: $tailCount sets found.\n";
open(VSS,">$dataD/vss") || die "could not open $dataD/vss";
my %seen;
my $procCount = 0;
foreach my $tail (keys(%tails))
{
    my $hits = $tails{$tail};
    if (@$hits <= 1)
    {
        $stats->Add(singletonTail => 1);
    } else {
        $_ = @$hits; print STDERR "sz hit = $_\n";
        my $longest = &longest($hits,\%seqs);
        print STDERR "Longest is $longest.\n";
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
                $stats->Add(vssOut => 1);
            } else {
                $stats->Add(vssSingleton => 1);
            }
        }
    }
    $procCount++;
    print STDERR "$procCount of $tailCount tail sets processed.\n";
}
close(VSS);
print STDERR "All done.\n" . $stats->Show();

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
