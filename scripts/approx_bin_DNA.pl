use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use gjoseqlib;

=head1 Use Genome Signatures to Bin DNA

    approx_bin_DNA Data < input-DNA

The Data directory must contain

    complete.genomes a 2-column table [GenomeId,GenomeName]
                     that defines the set of genomes for 
                     which signatures have been computed

    genome.signatures  a 2-column table  [Signature,GenomeId]
                     in Data.  Think of the GenomeIds as defining a set of
                     representative genomes, and the signatures as 8-mers that occur in
                     only one of the representative genomes (i.e., in a peg sequence of
                     the representative genome).

=head2 Input

The input is a fasta sequence of DNA fragments (usually contigs or reads)

=head2 Ouput

Executing the command produces a directory in which
there is a file for each genome in Data/complete.genomes.
The file will be in fasta format and will contain
the pieces of DNA believed to be in the partition represented
the the representative genome.

=cut

(
 (my $dataD = shift @ARGV) &&
 (my $bins  = shift @ARGV)
)
    || die "usgae: approx_bin_DNA Data Bins";

my $K=8;
if (-s "$dataD/K") {
    $K= &SeedUtils::file_head("$dataD/K", 1);
    chomp $K;
}

my($g2s,$s2g,$g2name) = &load_rep_sigs($dataD);

my $num_sigs_for_g = {};
foreach my $g (keys % { $g2s })
{
    my $sigs = $g2s->{$g};
    foreach my $sig (keys %$sigs)
    {
	++$num_sigs_for_g->{$g};
    }
}
open(BINS,">$bins") || die "could not open $bins";
while (my $tuple = &gjoseqlib::read_next_fasta_seq(\*STDIN))
{
    &process_tuple($tuple,$g2s,$s2g,$K,$num_sigs_for_g,\*BINS);
}
close(BINS);

sub process_tuple {
    my($tuple,$g2s,$s2g,$K,$num_sigs_for_g,$bins) = @_;
    
    my($id,$comment,$seq) = @$tuple;

    my $hits = {};
    &process_strand(\$seq,$g2s,$s2g,$K,$hits);
    my $rseq = &rev_comp($seq);
    &process_strand(\$rseq,$g2s,$s2g,$K,$hits);

    &process_hits($g2s,$s2g,$g2name,$id,$seq,$bins,$hits)
}

sub process_strand {
    my($seqP,$g2s,$s2g,$K,$hits) = @_;

    my $offset;
    for ($offset=0; ($offset <= 2); $offset++)
    {
	my $tmp = substr($$seqP,$offset);
	my $aa_seq = &SeedUtils::translate($tmp,undef,0);
	
	my $i;
	my $last = (length($aa_seq)- 1) - $K;
	for ($i=0; ($i <= $last); $i++)
	{
	    my $kmer = uc substr($aa_seq,$i,$K);
	    if (my $g = $s2g->{$kmer})
	    {
		$hits->{$g}->{$kmer} = 1;
	    }
	}
    }
}


sub process_hits {
    my($g2s,$s2g,$g2name,$id,$dna,$bins,$hits) = @_;
    
    my %num_hits_for_g;
    foreach my $g (keys % { $hits })
    {
	my $sigs = $hits->{$g};
	foreach my $sig (keys %$sigs)
	{
	    ++$num_hits_for_g{$g};
	}
    }
    
    my %frac_hits;
    foreach my $g (keys % { $g2s })
    {
	$frac_hits{$g} = $num_hits_for_g{$g} / $num_sigs_for_g{$g}
    }
    my @poss = sort { $frac_hits{$b} <=> $frac_hits{$a} } keys %frac_hits;
    my $g = $poss[0];
    next if ($frac_hits{$g} == 0.0);
    print $bins (join("\t", ($g,$id,$dna)), "\n");
    }
}


sub load_rep_sigs {
    my($dataD) = @_;

    my $g2s = {};
    my $s2g = {};
    open(SIGS,"<$dataD/genome.signatures")
	|| die "could not open $dataD/genome.signatures";
    while (defined($_ = <SIGS>))
    {
	if ($_ =~ /^(\S+)\t(\S+)/)
	{
	    my($sig,$gid) = ($1,$2);
	    $g2s->{$gid}->{$sig}++;
	    $s2g->{$sig} = $gid;
	}
    }
    close(SIGS);
    
    my $g2name;
    %$g2name = map { m/^(\S+)\t(\S.*\S)/ ? ($1 => $2) : () } &SeedUtils::file_read("$dataD/complete.genomes");
    
    return ($g2s,$s2g,$g2name);
}
