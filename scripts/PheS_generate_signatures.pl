use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use gjoseqlib;

=head1 Generate Signatures for PheS amino acid strings

    PheS_generate_signatures Data.PheS

The Data.PheS directory must contain

    complete.genomes a 2-column table [GenomeId,GenomeName]
                     that defines the set of genomes for 
                     which PheS signatures are computed

    6.1.1.20.fasta   which is a collection of PheS aa sequences.  These need 
                     to include all genomes in complete.genomes

=head2 Ouput

Executing the command creates a 2-column table [Signature,GenomeId]
in Data.PheS.  Think of the GenomeIds as defining a set of
representative genomes, and the signatures as 8-mers that occur in
only one of the representative genomes (i.e., in the PheS sequence of
the representative genome).

=cut

(my $dataD = shift @ARGV)
    || die "usage: pheS_generate_signatures Data";

my @tuples = gjoseqlib::read_fasta("$dataD/6.1.1.20.fasta");
(@tuples > 0) || die "Missing $dataD/6.1.1.20.fasta";
my %seqs = map { ($_->[1] => $_->[2]) } @tuples;
           
open(GENOMES,"<$dataD/complete.genomes")
    || die "Missing $dataD/complete.genomes";
my %g_names = map { ($_ =~ /^(\S+)\t(\S.*\S)$/) ? ($1 => $2) : () } <GENOMES>;
close(GENOMES);

my @gids    = sort { $a <=> $b } keys(%g_names);

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

# my $aa_alphabet = "ACDEFGHIKLMNPQRSTVWY";
my %kmers;
foreach my $gid (@gids)
{
    my $aa = uc $seqs{$gid};
    if (! $aa) { die "missing sequence for $gid" }
    my $i;
    my $last = (length($aa) - 1) - $K;
    for ($i = 0; ($i <= $last); $i++)
    {
	my $x = substr($aa,$i,$K);
	if ($x =~ /^[ACDEFGHIKLMNPQRSTVWY]{$K}$/o)
	{
	    $kmers{$x}->{$gid} = 1;
	}
    }
}

open(KMERS,">$dataD/PheS.signatures")
    || die "could not open $dataD/PheS.signatures";
foreach my $kmer (keys(%kmers))
{
    my $gidH = $kmers{$kmer};
    my @gids = keys(%$gidH);
    if (@gids == 1)
    {
	print KMERS join("\t",($kmer,$gids[0])),"\n";
    }
}
