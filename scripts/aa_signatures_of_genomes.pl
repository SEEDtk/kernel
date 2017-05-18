use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use ScriptUtils;
use Shrub;

=head1 Generate Protein Kmer Signatures for a Set of Genomes

    aa_signatures_of_genomes.pl [options] <genomes >signatures 2>log

This file generates signature protein kmers for genomes in a specified set. It is designed
to be applied to the output of L<PheS_evidence_of_reps>, but can be used for any set of
genomes by specifying the C<--filter=0> option.

=head2 Parameters

The standard input should be a tab-delimited file with genome IDs in the first column and scores in the second.
All genomes with scores above a certain level will be processed (default 60).

The command-line options are those found in L<ScriptUtils/ih_options> to specify the standard input, those
found in L<Shrub/script_options> to specify the database, plus the following.

=over 4

=item kmerSize

The kmer size for the signatures. The default is C<8>.

=item filter

The minimum score for which a genome should be considered representative. The default is C<60>. If C<0> is
specified, all genomes in the input file will be processed and no score column is necessary.

=back

=head2 Ouput

Executing the command creates a 2-column table [Signature,GenomeId]
on the standard output.  Think of the GenomeIds as defining a set of
representative genomes, and the signatures as 8-mers that occur in
only one of the representative genomes. Status messages will be
put on the standard error stream.

=cut

my $opt = ScriptUtils::Opts('dataDir', ScriptUtils::ih_options(), Shrub::script_options(),
        ['kmerSize|kmer|k=i', 'kmer size (default 8)', { default => 8 }],
        ['filter|f=i', 'minimum acceptable genome score', { default => 60 }]
);
my $K = $opt->kmersize;
# Read the genomes.
my $ih = ScriptUtils::IH($opt->input);
my @genomes;
my $filter = $opt->filter;
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^(\d+\.\d+)(?:\t(\d+))?/) { # Note that if the score is floating-point, the fractional part is ignored.
        my ($genome, $score) = ($1, $2 // 0); # Note that a missing or invalid score is turned into 0.
        if ($score >= $filter) {
            push @genomes, $genome;
        }
    }
}
print STDERR scalar(@genomes) . " input genomes found.\n";
my $shrub = Shrub->new_for_script($opt);
print STDERR "Connected to database.\n";
my $kmers = {};
foreach my $g (@genomes)
{
    print STDERR "Processing $g for $K-mers.\n";
    &load_kmers_for_genome($shrub, $g, $K, $kmers);
}

my ($kcount, $kept) = (0, 0);
foreach my $kmer (keys(%$kmers))
{
    my @genomes = keys(%{$kmers->{$kmer}});
    if (@genomes == 1)
    {
        print $kmer,"\t",$genomes[0],"\n";
        $kept++;
    }
    $kcount++;
    if ($kcount % 100000 == 0) {
        print STDERR "$kcount kmers processed. $kept kept.\n";
    }
}
print STDERR "$kcount total kmers. $kept kept.\n";

sub load_kmers_for_genome {
    my($shrub, $g,$K,$kmers) = @_;
    my $qh = $shrub->Get("Feature2Protein Protein",
                         'Feature2Protein(from-link) LIKE ?',
                         ["fig|$g.peg.%"], [qw(Protein(sequence))]);
    while (my $resultRow = $qh->Fetch()) {
        my $aa = $resultRow->PrimaryValue('Protein(sequence)');
        my $last = (length($aa) - 1) - $K;
        for (my $i = 0; ($i <= $last); $i++)
        {
            my $x = substr($aa,$i,$K);
            if ($x =~ /^[ACDEFGHIKLMNPQRSTVWY]{$K}$/o)
            {
                $kmers->{$x}->{$g} = 1;
            }
        }
    }
}