=head1 Build a Crude Representative Genome Set

    p3-crude-reps.pl [options] inDir outDir

This script will look at a directory built by L<p3-rep-prots.pl> and extract a list of representative genomes none of which have a similarity
score greater than the specified value.

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

For the input directory, the last two files are not used. This is considered a crude script because the C<rep_db.tbl> produced will connect
each genome to a representative, but not necessarily the closest, as a closer representative may appear later in the input.

=head2 Parameters

The positional parameters are the name of the input directory and the name of the output directory. The output directory will be turned into a full-blown
representative-genome directory with a table connecting each genome from the input directory to its representative.

The similarity score is the number of protein kmers in common.

The following command-line options are supported.

=over 4

=item kmer

The proposed kmer size for similarity comparison.  The default is C<8>.

=item clear

If specified, the output directory will be erased before any output is produced. The default is to leave existing files and overwrite any files already in place.

=item similarity

The maximum number of kmers in common for two genomes to be considered close. The default is C<100>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use FastA;
use Stats;
use Time::HiRes;
use Math::Round;

# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir',
        ["kmer|K|k=i", 'protein kmer size', { default => 8 }],
        ["clear", 'clear output directory'],
        ["similarity|sim|s=i", 'maximum similarity score', { default => 100 }],
        );
# Create the statistics object.
my $stats = Stats->new();
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
} elsif ($opt->clear) {
    print "Erasing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase $outDir: $!";
} else {
    print "Output directory is $outDir.\n";
}
# Create a new, blank representative-genome object.
my $repDb = RepGenomeDb->new(K => $opt->kmer, score => $opt->similarity);
print "Selected kmer size is " . $repDb->K . " for similarity " . $repDb->score . ".\n";
# Read all the input genomes.
print "Reading $inDir/complete.genomes.\n";
my %genomes;
open(my $kh, '<', "$inDir/complete.genomes") || die "Could not open genome file: $!";
while (! eof $kh) {
    my $line = <$kh>;
    if ($line =~ /^(\d+\.\d+)\s+(.+)/) {
        $genomes{$1} = $2;
        $stats->Add(genomeName => 1);
    }
}
close $kh; undef $kh;
print scalar(keys %genomes) . " found in file.\n";
# We will count the number of genomes processed and the number of representatives kept.
my ($gCount, $rCount) = (0, 0);
# Now we read in the FASTA file, keeping any that are not close to existing genomes.
print "Reading proteins.\n";
my $fh = FastA->new("$inDir/6.1.1.20.fasta");
my $start0 = time;
while ($fh->next) {
    # Get the genome ID.
    my $genome = $fh->id;
    if ($genome =~ /^fig\|(\d+\.\d+)/) {
        $genome = $1;
    }
    my $name = $genomes{$genome};
    if (! $name) {
        die "No name found for $genome.";
    }
    $stats->Add(proteinIn => 1);
    $gCount++;
    # Get the sequence.
    my $prot = $fh->left;
    # Find its representative (if any).
    my ($repID, $score) = $repDb->check_rep($prot);
    if ($repID) {
        # We are represented already.
        $repDb->Connect($repID, $genome, $score);
        $stats->Add(genomeConnected => 1);
    } else {
        # This is a new representative.
        $repDb->Add($genome, $name, $prot);
        $rCount++;
        print "New representative genome $genome: $name\n";
        $stats->Add(genomeChosen => 1);
    }
    if ($gCount % 5000 == 0) {
        print "$gCount genomes processed with $rCount representatives at " . Math::Round::nearest(0.01, (time - $start0) / $gCount) . " seconds/genome.\n";
    }
}
# All done. Write the output directory.
print "Writing output.\n";
$repDb->Save($outDir);
print "All done.\n" . $stats->Show();