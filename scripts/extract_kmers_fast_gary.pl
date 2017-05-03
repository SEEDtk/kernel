use strict;
use Data::Dumper;
use gjoseqlib;
use SeedUtils;
use ScriptUtils;
use Shrub;
use Shrub::Contigs;
use Time::HiRes;

my $opt = ScriptUtils::Opts('kmerSize', ScriptUtils::ih_options(), Shrub::script_options());
my $ih = ScriptUtils::IH($opt->input);
my $shrub = Shrub->new_for_script($opt);
my $k = shift @ARGV;
if (! $k) { $k=14; }
elsif ($k % 2) {
    die "K must be even";
}
my $genomes = 0;
my %all_kmers;
my $start0 = time;
my $starti = time;
while (defined($_ = <$ih>))
{
    if ($_ =~ /(\d+\.\d+)/)
    {
        my $g = $1;
        print STDERR "processing $g\n";
        my $contigs = Shrub::Contigs->new($shrub, $g);
        if (! $contigs) {
            print STDERR "Genome $g not found.\n";
        } else {
            my @tuples = $contigs->tuples;
            foreach my $tuple (@tuples)
            {
                my($gid,$comment,$seq) = @$tuple;

                my $s =&extract_kmers($seq,$k);
                foreach my $kmer (@$s) {
                    my $rev = SeedUtils::rev_comp($kmer);
                    if ($rev lt $kmer) {
                        $kmer = $rev;
                    }
                    my $x = $all_kmers{$kmer};
                    if (! $x)    {
                        $all_kmers{$kmer} = $g
                    } elsif ($x ne $g) {
                        $all_kmers{$kmer} = '*';
                    }
                }



#                my $n = length($seq);
#                for (my $i=0; ($i < ($n - $k)); $i++)
#                {
#                    my $kmer = lc substr($seq,$i,$k);
#                    if ($kmer =~ /^[acgt]*$/)
#                    {
#                        my $rev = SeedUtils::rev_comp($kmer);
#                        if ($rev lt $kmer) {
#                            $kmer = $rev;
#                        }
#                        my $x = $all_kmers{$kmer};
#                        if (! $x)    {
#                            $all_kmers{$kmer} = $g
#                        } elsif ($x ne $g) {
#                            $all_kmers{$kmer} = '*';
#                        }
#                    }
#                }

            }
            $genomes++;
            if ($genomes % 10 == 0) {
                my $duration = (time - $starti);
                $starti = time;
                print STDERR "$duration seconds per 10 genomes.\n";
            }
        }
    }
}
my $duration = (time - $start0);
print STDERR "$duration seconds per $genomes genomes.\n";
my ($found, $kept) = (0, 0);
while (my($kmer,$g) = each %all_kmers)
{
    $found++;
    if ($g ne '*')
    {
        print join("\t",($kmer,$g)),"\n";
        $kept++;
    }
}
$duration = (time - $start0);
print STDERR "$duration seconds total.\n";
print STDERR "$found kmers, $kept kept.\n";


sub extract_kmers {
    my($seq,$K) = @_;

    my $triples = int($K/2);
    my @kmers;
    my $pos = 0;
    my $last = (length($seq) - (3 * $triples));

    while ($pos <= $last)
    {
	push(@kmers, &extract_kmer(\$seq,$pos,$triples));
	$pos++;
    }
    return \@kmers;
}

#
#  $seq is a sequence (a contig)
#  $pos is an index of where to start pulling a k-mer
#  $triples is the number of "2of3" characters we pull.
#
#  The returned string is composed of $trples 2of3 chunks
#  Thus, a $triples value of 8 would pull 16 characters.
#
sub extract_kmer {
    my($seqP,$pos,$triples) = @_;

    my(@chars);
    my $i;
    for ($i=0; ($i < $triples); $i++)
    {
        my $kmer = lc substr($$seqP,$pos+($i*3),2);
        if ($kmer =~ /^[acgt]*$/) {
            push(@chars, $kmer);
            #push(@chars,lc substr($$seqP,$pos+($i*3),2));
        }
    }
    return join("",@chars);
}
