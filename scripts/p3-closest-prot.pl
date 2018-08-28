=head1 Find Closest protein

    p3-closest-prot.pl [options] fastaFile protein

Given a protein FASTA file, find the closest one to the specified input protein using kmer similarity. The
input protein can be expressed as a FASTA file name, a PATRIC feature ID, or an amino acid sequence string.
The output will be a single line containing the ID of the closest protein and the similarity score.

=head2 Parameters

The positional parameters are the name of a protein FASTA file and the input protein itself. The input protein can
be a FASTA file name (in which case the first protein will be used), a PATRIC feature ID, or a raw protein sequence.

The command-line options are as follows.

=over 4

=item kmerSize

The number of amino acids to use in the kmers generated for comparison purposes. The default is C<8>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenome;
use FastA;

# Get the command-line options.
my $opt = P3Utils::script_opts('fastaFile protein',
        ['kmerSize|kmersize|K=i', 'protein kmer size'],
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the positional parameters.
my ($fastaFile, $protein) = @ARGV;
if (! $fastaFile) {
    die "No protein FASTA file specified.";
} elsif (! -s $fastaFile) {
    die "Protein FASTA file $fastaFile missing or invalid.";
} elsif (! $protein) {
    die "No input protein specified.";
} elsif (-s $protein) {
    # Here the protein is probably a FASTA file.
    my $fh = FastA->new($fastaFile);
    if (! $fh->next) {
        die "$fastaFile appears to be empty.";
    } else {
        $protein = $fh->left;
    }
} elsif ($protein =~ /^fig\|/) {
    my ($protRecord) = $p3->query('genome_feature', ['select', 'aa_sequence'], ['eq', 'patric_id', $protein]);
    if (! $protRecord) {
        die "$protein not found in PATRIC.";
    } else {
        $protein = $protRecord->{aa_sequence};
    }
}
# Create the kmer hash for the input protein.
my $rep = RepGenome->new('input', prot => $protein, K => $opt->kmersize);
# Loop through the protein FASTA, remembering the closest match.
my $best = '<none>';
my $score = 0;
my $fh = FastA->new($fastaFile);
while ($fh->next) {
    my $newProt = $fh->left;
    my $newScore = $rep->check_genome($newProt);
    if ($newScore > $score) {
        $best = $fh->id;
        $score = $newScore;
    }
}
# Write what we found.
print "$best\t$score\n";
