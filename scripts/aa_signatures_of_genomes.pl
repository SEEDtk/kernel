use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use gjoseqlib;

=head1 Generate AA Signatures for a Set of Genomes

    aa_signatures_of_genomes.pl Data

The Data directory must contain

    complete.genomes a 2-column table [GenomeId,GenomeName]
                     that defines the set of genomes for
                     which aa signatures are computed

    K                typically 8-mers are used

=head2 Ouput

Executing the command creates a 2-column table [Signature,GenomeId]
in Data.  Think of the GenomeIds as defining a set of
representative genomes, and the signatures as 8-mers that occur in
only one of the representative genomes

=cut

(my $dataD = shift @ARGV)
    || die "usage: aa_signatures_of_genomes Data";

my $K = 8;
if (open(K,"<$dataD/K"))
{
    $_ = <K>;
    close(K);
    if ($_ =~ /^(\d+)/)
    {
        $K = $1;
    }
}

open(GENOMES,"<$dataD/complete.genomes")
    || die "could not open $dataD/complete.genomes";
my @genomes = map { ($_ =~ /^(\S+)/) ? $1 : () } <GENOMES>;
close(GENOMES);

open(SIG,">$dataD/genome.signatures")
    || die "could not open $dataD/genome.signatures";

my $kmers = {};
foreach my $g (@genomes)
{
    &load_kmers_for_genome($g,$K,$kmers);
}
foreach my $kmer (keys(%$kmers))
{
    my @genomes = keys(%{$kmers->{$kmer}});
    if (@genomes == 1)
    {
        print SIG $kmer,"\t",$genomes[0],"\n";
    }
}
close(SIG);

sub load_kmers_for_genome {
    my($g,$K,$kmers) = @_;

    my $peg_sequences = &get_peg_sequences_for_genome($g);

    # my $aa_alphabet = "ACDEFGHIKLMNPQRSTVWY";

    foreach my $peg (keys(%$peg_sequences))
    {
        my $aa = uc $peg_sequences->{$peg};
        if (! $aa) { die "missing sequence for $peg" }
        my $i;
        my $last = (length($aa) - 1) - $K;
        for ($i = 0; ($i <= $last); $i++)
        {
            my $x = substr($aa,$i,$K);
            if ($x =~ /^[ACDEFGHIKLMNPQRSTVWY]{$K}$/o)
            {
                $kmers->{$x}->{$g} = 1;
            }
        }
    }
}