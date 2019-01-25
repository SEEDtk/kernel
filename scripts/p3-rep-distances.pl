=head1 Compute Representative-Genome Distances

    p3-rep-distances.pl [options] repDir

This script will read a FASTA file of seed protein sequences and return the distance of each to all representative genomes
represented by it.  So, if a particular genome is sufficiently close to three genomes in the representative set, three
distances will be returned.  This will be done for each genome in the input FASTA, which can be a very lengthy process.

=head2 Parameters

The positional parameter is the name of the representative-genome directory as described in L<RepGenomeDb>.

The standard input should contain a FASTA file.  The ID field should be a genome or feature ID and the sequence field the
amino acid sequence of the same seed protein used to build the representative-genome directory.  Feature IDs will be
converted automatically to genome IDs.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The standard output will contain a tab-delimited file with three columns-- (0) first genome ID, (1) second genome ID, and
(2) distance from 0 to 1.

The additional command-line options are as follows.

=over 4

=item verbose

Status messages will be written to the standard error output.

=item comprehensive

Distances between all the pairs of representative genomes will also be output.  NOTE: to get just the matrix, use this option
with an empty input file.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use RepGenome;
use FastA;
use Stats;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir', P3Utils::ih_options(),
        ['verbose|debug|v', 'write status messages to STDERR'],
        ['comprehensive|matrix|X', 'also output full representative-genome distances']
        );
my $debug = $opt->verbose;
my $matrix = $opt->comprehensive;
# Initialize our statistical tracking.
my $stats = Stats->new(qw(genomeIn distanceOut));
# Get the input directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Invalid or missing input directory $repDir.";
}
# Load the database.
print STDERR "Loading database from $repDir.\n" if $debug;
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Get the kmer size and minimum similarity scores for this database.
my $K = $repDB->K;
my $minScore = $repDB->score;
# Open the input file.
my $ih = P3Utils::ih($opt);
my $fastaH = FastA->new($ih);
# This will count genomes in for progress messages.
my $count = 0;
# Write the output header.
P3Utils::print_cols(['genome1', 'genome2', 'distance']);
# Now we loop through the input FASTA file.
while ($fastaH->next()) {
    my ($genome) = ($fastaH->id =~ /(\d+\.\d+)/);
    die "Invalid sequence ID " . $fastaH->id . " in input file." if ! $genome;
    # Is this genome in the representative set?
    if ($repDB->is_rep($genome)) {
        # Yes.  Skip it.
        $stats->Add(repFound => 1);
    } else {
        # No, we must process this genome.  We start with its sequence, then get the nearby representatives.
        my $protSeq = $fastaH->left;
        my $repHash = $repDB->list_reps($protSeq, $minScore);
        my @others = keys %$repHash;
        if (scalar(@others) > 1) {
            $stats->Add(multiRepGenome => 1);
        } elsif (! @others) {
            $stats->Add(unreppedGenome => 1);
            print STDERR "WARHING: Genome $genome is unrepresented.\n" if $debug;
        }
        # Now we have the representative genomes near ours.  Create a repGenome object for our new genome.
        my $repObject = RepGenome->new($genome, prot => $protSeq, K => $K);
        # Use it to compute the distances.
        for my $other (@others) {
            my $otherObject = $repDB->rep_object($other);
            my $distance = $repObject->distance($otherObject);
            P3Utils::print_cols([$genome, $other, $distance]);
            $stats->Add(distanceOut => 1);
        }
    }
    $count++;
    print STDERR "$count input genomes processed.\n" if $debug && $count % 100 == 0;
}
print STDERR "$count total genomes read from input.\n" if $debug;
# Check for the comprehensive option.
if ($matrix) {
    # Here the user wants the matrix of distances between representatives.
    my $genomeList = $repDB->rep_list();
    my $total = scalar @$genomeList;
    print STDERR "Creating representative distance matrix.  $total genomes in list.\n" if $debug;
    $count = 0;
    # Start at the beginning of the list.
    my $genome = shift @$genomeList;
    while (scalar @$genomeList) {
        my $repObject = $repDB->rep_object($genome);
        $count++;
        $stats->Add(repGenomeIn => 1);
        print STDERR "Processing $genome ($count of $total).\n" if $debug;
        for my $other (@$genomeList) {
            my $otherObject = $repDB->rep_object($other);
            my $distance = $repObject->distance($otherObject);
            P3Utils::print_cols([$genome, $other, $distance]);
            $stats->Add(repDistanceOut => 1);
        }
        # Get the next genome.
        $genome = shift @$genomeList;
    }
}
print STDERR "All done.\n" . $stats->Show();
