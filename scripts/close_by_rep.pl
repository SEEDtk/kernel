use gjoseqlib;
use strict;
use Data::Dumper;

### constants
my $K     = 8;
my $min   = 100;
my $Nbest = 1;
#
# $K will probably always be 8, but who knows.  An
# analysis would be a good idea.
#
# $min is the minimum number of kmers in common for us
#      consider a sequence as close enough
#
# $Nbest is the number of close genome Ids you want to get back
#      Now the output is Nbest lines, each containing a score and an id (sorted)
###
# my $usage = "usage: close_by_rep PheSs SinglePheS [Nbest]";
#
#
# PheSs is a fasta file of PheS sequences.  We seek the closest Nbest
# of these to the single sequence in the file given by
# the second argument

my $PheS_file;
my $org_seq_file;
my $usage;
(
 ($PheS_file    = shift @ARGV) &&
 ($org_seq_file = shift @ARGV)
)
    || die $usage;

$Nbest = shift @ARGV; if (! $Nbest ) { $Nbest = 1 }

# The rep_seqs (from the file given as the first argument,
# will usually contain a set of representative genes,
# but not always.  They are just a set of PheS sequences,
# and you want to find the closest.

my @rep_seqs = &gjoseqlib::read_fasta($PheS_file);
(@rep_seqs > 1) || die "bad file of reps";

####################
# This second argument should be a fasta file containing just
# one sequence, the query sequence.
my @tmp     = &gjoseqlib::read_fasta($org_seq_file);
my $seq1    = $tmp[0]->[2];
((@tmp == 1) && ($seq1 && (length($seq1) > $K)))
    || die "bad file containing singleton";
####################
#
# We begin by computing a hash of the Kmers that occur
# in the query sequence

my %seq1H;
my $i1;

for ($i1=0; ($i1 < (length($seq1) - $K)); $i1++)
{
    my $kmer1 = uc substr($seq1,$i1,$K);
    $seq1H{$kmer1} = 1;
}
#### We now have $seq1H as a hash of kmers from $seq1
####################

# Now we go through the rep_seqs (remember, these could be
# any set of genomes)
#
# for each rep genome, we just count the hits against the query
# kmers. We keep a running track of the rep_seq that has
# a maximum overlap.
# When we get through the pass thru the rep_seqs.
# we just print the one with the maximum overlap
# against Kmers of the query sequence.

# During the pass through the reps we maintain
#
#     best        the best genome candidate
#
#     best_sofar The overlap against the best candidate

my @best_sofar = ();

foreach my $tuple (@rep_seqs)
{
    my($id,$comment,$seq2) = @$tuple;
    my %hits;
    my $i2;
    my $max = length($seq2) - $K;

    for ($i2=0; ($i2 < $max); $i2++)
    {
        my $kmer2 = uc substr($seq2,$i2,$K);
        if ($seq1H{$kmer2})
        {
            $hits{$kmer2} = 1;
        }
    }

    my $hitsN = scalar keys(%hits);
    if (($hitsN && ($hitsN >= $min)) &&
        ((@best_sofar == 0) ||
         ($hitsN > $best_sofar[-1]->[0])))
    {
        my $i3;
        for ($i3=0; ($hitsN < $best_sofar[$i3]->[0]); $i3++) {}
        my $remove = (@best_sofar == $Nbest) ? 1 : 0;
        splice(@best_sofar,$i3,0,[$hitsN,$id]);
        pop @best_sofar if $remove;
    }
}

foreach my $tuple (@best_sofar)
{
    print join("\t",@$tuple),"\n";
}


