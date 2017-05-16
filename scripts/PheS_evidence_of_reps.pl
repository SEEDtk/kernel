use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use gjoseqlib;

=head1 Use PheS Signatures to Predict Presence of Representative Genomes

    PheS_evidence_of_reps Data.PheS < input-DNA

The Data.PheS directory must contain

    complete.genomes a 2-column table [GenomeId,GenomeName]
                     that defines the set of genomes for 
                     which PheS signatures are computed

    6.1.1.20.fasta   which is a collection of PheS aa sequences.  These need 
                     to include all genomes in complete.genomes

    PheS.signatures  a 2-column table  [Signature,GenomeId]
                     in Data.PheS.  Think of the GenomeIds as defining a set of
                     representative genomes, and the signatures as 8-mers that occur in
                     only one of the representative genomes (i.e., in the PheS sequence of
                     the representative genome).

=head2 Input

The input is a fasta sequence of DNA fragments (usually contigs or reads)

=head2 Ouput

Executing the command produces a report suggesting which
of the representative genomes are present in the input file

=cut

(my $dataD = shift @ARGV)
    || die "usgae: pheS_evidence_of_reps Data";

my $K=8;
if (-s "$dataD/K") {
    $K= &SeedUtils::file_head("$dataD/K", 1);
    chomp $K;
}

my($g2s,$s2g,$g2name) = &load_PheS_sigs($dataD);

my $hits = {};
while (my $tuple = &gjoseqlib::read_next_fasta_seq(\*STDIN))
{
    my $seq = $tuple->[2];
    &process_seq(\$seq,$g2s,$s2g,$K,$hits)
}

&process_hits($g2s, $s2g, $g2name, $hits);

sub process_hits {
    my($g2s, $s2g, $g2name, $hits) = @_;
    
    my %num_sigs_for_g;
    foreach my $g (keys % { $g2s })
    {
	my $sigs = $g2s->{$g};
	foreach my $sig (keys %$sigs)
	{
	    ++$num_sigs_for_g{$g};
	}
    }
    
#     print STDOUT map { (join("\t", @$_), "\n")
#     } sort { $b->[1] <=> $a->[1]
#     } map { [$_, (scalar keys % { $hits->{$_} })] } keys % { $hits };
#     die "aborted";
    
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
    
    foreach my $g (sort { $frac_hits{$b} <=> $frac_hits{$a} } keys %frac_hits)
    {
	next if ($frac_hits{$g} == 0.0);
	print STDOUT (join("\t", ($g, 100.0*$frac_hits{$g}, $g2name->{$g})), "\n");
    }
}


sub process_seq {
    my($seqP,$g2s,$s2g,$K,$hits) = @_;
    
    &process_strand($seqP,$g2s,$s2g,$K,$hits);
    &rev_comp($seqP);    #...rev_comp value is returned "In-place"...
    &process_strand($seqP,$g2s,$s2g,$K,$hits);
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

sub load_PheS_sigs {
    my($dataD) = @_;

    my $g2s = {};
    my $s2g = {};
    open(SIGS,"<$dataD/PheS.signatures")
	|| die "could not open $dataD/PheS.signatures";
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
