
=head1 Count Representative Genomes

    p3-count-reps.pl [options] inDir score

This script looks at an input list of PATRIC genome IDs and determines the number of representative genomes for each one
at the specified distance or closer.

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

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain PATRIC genome IDs in the key column.

Additional command-line options are those given in L<P3Utils/col_options> (to select the input column).

=head2 Output Files

The standard output will be tab-delimited, with each record containing (0) a genome ID, (1) the genome name, and (2) the number of representatives
found at the given distance.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir score', P3Utils::col_options(), P3Utils::ih_options(),
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
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
# Create the PATRIC filter and column clauses for genome queries.
my @filter = (['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']);
my @cols = qw(genome_id genome_name aa_sequence);
# Create the database from the input directory.
my $repDB = RepGenomeDb->new_from_dir($inDir);
# Save the parameters.
my $K = $repDB->K();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Print the output headers.
P3Utils::print_cols([qw(genome_id genome_name count)]);
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    my @genomes = map { $_->[0] } @$couplets;
    # Get the name and seed protein for each genome.
    my $resultList = P3Utils::get_data_keyed($p3, 'feature', \@filter, \@cols, \@genomes, 'genome_id');
    # The resultList entries are in the form [$genome, $name, $prot]. Get the longest
    # protein for each genome.
    my (%results);
    for my $result (@$resultList) {
        my ($genome, $name, $prot) = @$result;
        # Check the protein.
        if ($prot) {
            # Add the protein length to the result array.
            my $protLen = length $prot;
            push @$result, $protLen;
            if (! exists $results{$genome}) {
                # It's the first protein, so store it.
                $results{$genome} = $result;
            } elsif ($protLen > $results{$genome}[3]) {
                # It's a better protein, so keep it.
                $results{$genome} = $result;
            }
        }
    }
    # Now loop through the genomes, checking for representatives.
    for my $genome (@genomes) {
        my $genomeData = $results{$genome};
        if ($genomeData) {
            my (undef, $name, $prot) = @$genomeData;
            my $count = $repDB->count_rep($prot, $score);
            P3Utils::print_cols([$genome, $name, $count]);
        }
    }
}

