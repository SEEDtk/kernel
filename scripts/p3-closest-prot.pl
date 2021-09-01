=head1 Find Closest protein

    p3-closest-prot.pl [options] fastaFile protein

Given a protein FASTA file, find the closest one to the specified input protein using kmer similarity. The
input protein can be expressed as a FASTA file name, a BV-BRC feature ID, or an amino acid sequence string.
The output will be a single line containing the ID of the closest protein and the similarity score.

=head2 Parameters

The positional parameters are the name of a protein FASTA file and the input protein itself. The input protein can
be a FASTA file name (in which case the first protein will be used), a BV-BRC feature ID, or a raw protein sequence.

The command-line options are as follows.

=over 4

=item kmerSize

The number of amino acids to use in the kmers generated for comparison purposes. The default is C<8>.

=item count

The number of matches to return. The default is C<1>, returning a single match.

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
        ['count|N=i', 'number of matches to return', { default => 1 }],
        );
# Get access to BV-BRC.
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
    my $fh = FastA->new($protein);
    if (! $fh->next) {
        die "$fastaFile appears to be empty.";
    } else {
        $protein = $fh->left;
    }
} elsif ($protein =~ /^fig\|/) {
    my ($protRecord) = P3Utils::get_data($p3, feature => [['eq', 'patric_id', $protein]], ['patric_id', 'aa_sequence']);
    if (! $protRecord) {
        die "$protein not found in BV-BRC.";
    } else {
        $protein = $protRecord->[1];
    }
}
# Create the kmer hash for the input protein.
my $rep = RepGenome->new('input', prot => $protein, K => $opt->kmersize);
# The following list will contain [id, score] tuples, sorted from best score to worst.
my @best;
# This is the number of proteins we want.
my $count = $opt->count;
# Loop through the protein FASTA, remembering the closest matches.
my $fh = FastA->new($fastaFile);
while ($fh->next) {
    my $prot = $fh->left;
    my $score = $rep->check_genome($prot);
    if ($score) {
        # Find where this score goes in the list.
        my $i;
        for ($i = $#best; $i >= 0 && $score > $best[$i][1]; $i--) {}
        $i++;
        if ($i < $count) {
            splice @best, $i, 0, [$fh->id, $score];
            if (scalar(@best) > $count) { pop @best; }
        }
    }
}
# Write what we found.
print "id\tscore\n";
print map { "$_->[0]\t$_->[1]\n" } @best;

