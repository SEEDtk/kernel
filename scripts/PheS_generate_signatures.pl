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

$| = 1;
(my $dataD = shift @ARGV)
    || die "usage: pheS_generate_signatures Data";
print "Reading FASTA file.\n";
my @tuples = gjoseqlib::read_fasta("$dataD/6.1.1.20.fasta");
(@tuples > 0) || die "Missing $dataD/6.1.1.20.fasta";
my %seqs;
for my $tuple (@tuples) {
    my ($id,$comment,$aa) = @$tuple;
    if ($id =~ /(\d+\.\d+)/) {
        $seqs{$1} = uc $aa;
    } elsif ($comment =~ /^(\d+\.\d+)/) {
        $seqs{$1} = uc $aa;
    }
}
print "Reading genome list.\n";
open(GENOMES,"<$dataD/complete.genomes")
    || die "Missing $dataD/complete.genomes";
my %g_names = map { ($_ =~ /^(\S+)\t(\S.*\S)$/) ? ($1 => $2) : () } <GENOMES>;
close(GENOMES);

my @genomes    = sort { $a <=> $b } keys(%g_names);
my $total = scalar @genomes;
print "$total genomes found.\n";
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
print "Kmer size is $K.\n";
# my $aa_alphabet = "ACDEFGHIKLMNPQRSTVWY";
my %kmers;
my $count = 0;
foreach my $gid (@genomes)
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
    $count++;
    print "$count of $total genomes processed.\n" if $count % 100 == 0;
}

open(KMERS,">$dataD/PheS.signatures")
    || die "could not open $dataD/PheS.signatures";
my $kTotal = scalar keys %kmers;
my $found = 0;
$count = 0;
my %sigs;
foreach my $kmer (keys(%kmers))
{
    my $gidH = $kmers{$kmer};
    my @gids = keys(%$gidH);
    if (@gids == 1)
    {
        print KMERS join("\t",($kmer,$gids[0])),"\n";
        $found++;
        $sigs{$gids[0]}++;
    }
    $count++;
    print "$count of $kTotal kmers processed. $found signatures found.\n" if $count % 1000 == 0;
}
print "All done. $found signatures found in $total genomes.\n\n";
print "genome_id\tcount\tname\n";
for my $gid (@genomes) {
    my $count = $sigs{$gid} // 0;
    my $name = $g_names{$gid} // '<unknown>';
    print "$gid\t$count\t$name\n";
}