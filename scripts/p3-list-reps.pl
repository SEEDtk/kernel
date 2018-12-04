=head1 List Representative Genomes For a Protein

    p3-list-reps.pl [options] inDir score

This script takes as input a FASTA file and lists the representative genomes at the specified distance or closer, along with
the associated score.

The script operates on a representative server directory. The directory contains three files of paramount interest.

=over 4

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

=head2 Parameters

The positional parameters are the input directory and the minimum similarity score. The input directory must contain the above three files.

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain the FASTA file with seed proteins (Phenylalanyl
tRNA-synthetase protein alpha chain) for multiple genomes.


=head2 Output Files

The standard output will be tab-delimited, with each record containing (0) a genome ID, (1) the genome name, and (2) the similarity score.

=cut

use strict;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use FastA;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir score', P3Utils::ih_options(),
        );
# Verify the parameters.
my ($inDir, $score) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! -f "$inDir/6.1.1.20.fasta" || ! -f "$inDir/complete.genomes") {
    die "$inDir does not appear to contain a representative-genomes database.";
}
if (! $score) {
    die "No minimum score specified.";
} elsif ($score =~ /\D/) {
    die "Invalid score $score.";
}
# Create the database from the input directory.
my $repDB = RepGenomeDb->new_from_dir($inDir);
# Save the parameters.
my $K = $repDB->K();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Print the output headers.
P3Utils::print_cols([qw(id rep_id genome_name similarity)]);
# The low-scoring genome output lines are kept in here.
my @outliers;
# Read the proteins.
my $reader = FastA->new($ih);
while ($reader->next) {
    my $prot = $reader->left;
    my ($repID, $score1) = $repDB->find_rep($prot);
    my $repName = '';
    if ($repID) {
        $repName = $repDB->rep_object($repID)->name;
    } else {
        $repID = '<none>';
    }
    my $output = [$reader->id, $repID, $repName, $score1];
    if ($score1 >= $score) {
        P3Utils::print_cols($output);
    } else {
        push @outliers, $output;
    }
}
# Write the outliers, if any.
if (@outliers) {
    print "\n";
    for my $outlier (@outliers) {
        P3Utils::print_cols($outlier);
    }
}