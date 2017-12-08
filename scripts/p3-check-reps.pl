
=head1 Check For Representative Genomes

    p3-check-reps.pl [options] inDir outDir

This script looks at an input list of PATRIC genome IDs and determines the representative genome for each one.

The script operates on a representative server directory. The directory contains four files of paramount interest.

=over 4

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item rep_db.tbl

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID, (1) the ID of the genome's
representative, and (2) the similarity number.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

=head2 Parameters

The positional parameters are the input directory and the output directory. The input directory must contain the above three files.

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain PATRIC genome IDs in the key column.

Additional command-line options are those given in L<P3Utils/col_options> (to select the input column).

=head2 Output Files

The standard output will contain a progress report.

The output directory will contain an updated version of the input directory with all the genomes represented.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir', P3Utils::col_options(), P3Utils::ih_options(),
        );
# Create the statistics object.
my $stats = Stats->new();
# Get access to PATRIC.
print "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# Verify the parameters.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! -f "$inDir/6.1.1.20.fasta" || ! -f "$inDir/complete.genomes") {
    die "$inDir does not appear to contain a representative-genomes database.";
} else {
    print "Input database in $inDir.\n";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} else {
    print "Output directory is $outDir.\n";
}
# Create the PATRIC filter and column clauses for genome queries.
my @filter = (['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']);
my @cols = qw(genome_id genome_name aa_sequence aa_len);
# Create the database from the input directory.
print "Creating database from $inDir.\n";
my $repDB = RepGenomeDb->new_from_dir($inDir, verbose => 1);
# Save the parameters.
my $K = $repDB->K();
my $minScore = $repDB->score();
print "Kmer size is $K and minimum similarity is $minScore.\n";
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# This will be a list of RepGenome objects for the un-represented genomes.
my @outliers;
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    my @genomes = map { $_->[0] } @$couplets;
    # We need to isolate the genomes for which we don't already have a representative.
    my @lost;
    for my $genome (@genomes) {
        $stats->Add(genomeIn => 1);
        my ($repID, $score) = $repDB->check_rep($genome);
        if (! $repID) {
            push @lost, $genome;
            $stats->Add(genomeNotFound => 1);
        } else {
            $stats->Add(genomeFound => 1);
        }
    }
    if (@lost) {
        print scalar(@lost) . " genomes in this batch are not yet in the database.\n";
        # Ask PATRIC for the names and identifying proteins of the un-represented genomes.
        my $resultList = P3Utils::get_data_keyed($p3, 'feature', \@filter, \@cols, \@lost, 'genome_id');
        # The resultList entries are in the form [$genome, $name, $prot, $protLen]. Get the longest
        # protein for each genome.
        my %results;
        for my $result (@$resultList) {
            my ($genome, $name, $prot, $protLen) = @$result;
            if (! exists $results{$genome}) {
                $results{$genome} = $result;
            } else {
                $stats->Add(redundantProt => 1);
                print "WARNING: $genome $name has a redundant identifying protein.\n";
                if ($protLen > $results{$genome}[3]) {
                    # It's a better protein, so keep it.
                    $stats->Add(redundantProt => 1);
                    $results{$genome} = $result;
                }
            }
        }
        # Now loop through the genomes, checking for representatives.
        for my $genome (keys %results) {
            print "Checking $genome.\n";
            my (undef, $name, $prot) = @{$results{$genome}};
            my ($repID, $score) = $repDB->find_rep($prot);
            if ($score >= $minScore) {
                print "$genome $name assigned to $repID with similarity $score.\n";
                $repDB->Connect($repID, $genome, $score);
                $stats->Add(genomeConnected => 1);
            } else {
                print "$genome $name is an outlier. Best score was $score.\n";
                my $repGenome = RepGenome->new($genome, name => $name, prot => $prot, K => $K);
                push @outliers, $repGenome;
                $stats->Add(genomeOutlier => 1);
            }
        }
    }
}
# Now we need to process the outliers.
my $nOutliers = scalar @outliers;
print "$nOutliers outlier genomes found.\n";
# Loop through the outliers, computing a cross-reference. The cross-reference is a two-dimensional table
# of similarity scores, indexed by the genome's position in the outlier list.
my $xref = [];
my ($i, $j);
for ($i = 0; $i < $nOutliers; $i++) {
    my $repGenome = $outliers[$i];
    print "Cross-referencing " . $repGenome->id() . "\n";
    # Copy the scores already computed.
    for ($j = 0; $j < $i; $j++) {
        $xref->[$i][$j] = $xref->[$j][$i];
    }
    # Compute the remaining scores.
    for ($j = $i; $j < $nOutliers; $j++) {
        my $prot = $outliers[$j]->prot();
        $xref->[$i][$j] = $repGenome->check_genome($prot);
    }
}
# Now we process the scores until we run out of outliers.
while (@outliers) {
    print "Searching for most popular genome among $nOutliers outliers.\n";
    # Find the genome with the most similarities.
    my ($b, $best) = (0, count($xref->[0], $minScore));
    for ($i = 1; $i < $nOutliers; $i++) {
        my $count = count($xref->[$i], $minScore);
        if ($count > $best) {
            $b = $i; $best = $count;
        }
    }
    my $bestGenome = $outliers[$b];
    my $bestID = $bestGenome->id();
    my $bestName = $bestGenome->name();
    print "Best genome is $bestID $bestName with $best neighbors.\n";
    # This list will contain the indices of the outliers we are keeping; that is, the
    # ones not similar to the best one.
    my @keepers;
    my $connected = 0;
    # Loop through the outliers, connecting the similar genomes to the best one.
    for ($i = 0; $i < $nOutliers; $i++) {
        my $score = $xref->[$b][$i];
        if ($score >= $minScore) {
            $bestGenome->AddGenome($outliers[$i]->id(), $score);
            $stats->Add(genomeConnected => 1);
            $connected++;
        } else {
            push @keepers, $i;
        }
    }
    print "$connected genomes represented by $bestID.\n";
    # Add the best genome to the database.
    $repDB->AddRepObject($bestGenome);
    $stats->Add(newRepAdded => 1);
    # Form the new outliers list.
    @outliers = map { $outliers[$_] } @keepers;
    $nOutliers = scalar @outliers;
    # Form the new xref. We map the old indices to the new ones, which eliminates the
    # genomes we connected above.
    my @newXref;
    print "Reorganizing the remaining $nOutliers genomes.\n";
    for ($i = 0; $i < $nOutliers; $i++) {
        for ($j = 0; $j < $nOutliers; $j++) {
            $newXref[$i][$j] = $xref->[$keepers[$i]][$keepers[$j]];
        }
    }
    $xref = \@newXref;
}
# All the outliers have been added to the database.
print "Saving to $outDir.\n";
$repDB->Save($outDir);
print "All done.\n" . $stats->Show();


# Count the good similarities in a list of scores.
sub count {
    my ($list, $minScore) = @_;
    my $retVal = 0;
    for my $score (@$list) {
        if ($score >= $minScore) {
            $retVal++;
        }
    }
    return $retVal;
}