use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use FastA;
use FastQ;
use ScriptUtils;

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

The input is a file of DNA fragments (usually contigs (FASTA) or reads (FASTQ)).

=head2 Parameters

The positional parameters should be the name of the working directory and optionally one or two input file names.
Use one file name for FASTA or interlaced FASTQ input. Use two file names for paired-end FASTQ input. If only one
input file is being used, you can omit all file names to use the standard input.

The following command-line parameters are supported.

=over 4

=item kmerSize

The kmer size for the input signatures. The default is C<8>.

=item fasta

If specified, the input is presumed to be in FASTA format instead of FASTQ format.

=back

=head2 Output

Executing the command produces a report suggesting which
of the representative genomes are present in the input file

=cut

my $opt = ScriptUtils::Opts('DataDir file1 file2',
        ['fasta', 'input is in FASTA format'],
        ['kmerSize|kmer|k=i', 'kmer size (default 8)', { default => 8 }]);
(my $dataD = shift @ARGV)
    || die "usgae: pheS_evidence_of_reps Data";
my ($file1, $file2) = @ARGV;
my $count = 0;
my $hitCount = 0;
my $K= $opt->kmersize;

my($g2s,$s2g,$g2name) = &load_PheS_sigs($dataD);

my $reader;
if ($opt->fasta) {
    $file1 //= \*STDIN;
    $reader = FastA->new($file1);
} else {
    $file1 //= \*STDIN;
    $reader = FastQ->new($file1, $file2);
}
my $hits = {};
while ($reader->next())
{
    my $seq = $reader->left;
    &process_seq(\$seq,$g2s,$s2g,$K,$hits);
    my $seq2 = $reader->right;
    &process_seq(\$seq,$g2s,$s2g,$K,$hits)

}
print STDERR "Analyzing the hits.\n";
&process_hits($g2s, $s2g, $g2name, $hits);

sub process_hits {
    my($g2s, $s2g, $g2name, $hits) = @_;

    my %num_sigs_for_g;
    foreach my $g (keys % { $g2s })
    {
        print STDERR "Counting hits for $g.\n";
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
        print STDERR "Summarizing hits for $g.\n";
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
        print STDOUT (join("\t", ($g, 100.0*$frac_hits{$g}, $num_hits_for_g{$g}, $g2name->{$g})), "\n");
    }
}


sub process_seq {
    my($seqP,$g2s,$s2g,$K,$hits) = @_;

    &process_strand($seqP,$g2s,$s2g,$K,$hits);
    &rev_comp($seqP);    #...rev_comp value is returned "In-place"...
    &process_strand($seqP,$g2s,$s2g,$K,$hits);
    $count++;
    if ($count % 100000 == 0) {
        print STDERR "$count sequences processed. $hitCount hits.\n";
    }
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
                $hitCount++;
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
