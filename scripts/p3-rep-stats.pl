=head1 Display Statistics on Representative Genomes

    p3-rep-stats.pl [options] repDir

This script produces a simple report on the number of genomes represented by each genome in a representative-genome database.

Status messages will be written to the standard error output.

=head2 Parameters

The single positional parameter is the name of the directory containing the database. The files in the directory are
described in L<RepGenomeDb>.

The command-line options are as follows.

=over 4

=item matrix

If specified, the name of a file to contain a report on distances between the representatives themselves.  The file
will contain three tab-delimited columns-- (0) the first genome ID, (1) the second genome ID, and (2) the distance.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir',
        ['matrix|X=s', 'output file for distance matrix']);
# Get the input directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Invalid or missing input directory $repDir.";
}
# Load the database.
print STDERR "Loading database from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Get the list of genomes.
print STDERR "Counting set sizes.\n";
my $genomeList = $repDB->rep_list();
# Create a hash mapping each genome to its number of representatives.
my %repHash;
for my $genomeID (@$genomeList) {
    my $repGenome = $repDB->rep_object($genomeID);
    my $repList = $repGenome->rep_list();
    $repHash{$genomeID} = scalar @$repList;
}
# Sort the hash by representative count.
my @genomeIDs = sort { $repHash{$b} <=> $repHash{$a} } @$genomeList;
# Write the report.
print STDERR "Writing report.\n";
print join("\t", qw(id count protLen name)) ."\n";
for my $genome (@genomeIDs) {
    my $repGenome = $repDB->rep_object($genome);
    my $protLen = length $repGenome->prot();
    my $name = $repGenome->name();
    print join("\t", $genome, $repHash{$genome}, $protLen, $name) . "\n";
}
# Check for the matrix.
my $outFile = $opt->matrix;
if ($outFile) {
    open(my $oh, '>', $outFile) || die "Could not open matrix file $outFile: $!";
    print $oh join("\t", qw(id1 id2 distance)) . "\n";
    # These will track our progress.
    my ($num, $total) = (1, scalar @genomeIDs);
    # These will be used for statistical metrics.
    my ($min, $max, $sum, $count) = (1, 0, 0, 0);
    # Get the first genome.
    my $genome = shift @genomeIDs;
    while (@genomeIDs) {
        # Get the first genome's object.
        my $repGenome = $repDB->rep_object($genome);
        # Compare it to the other genomes.
        print STDERR "Processing matrix row for $genome ($num of $total).\n";
        for my $genome2 (@genomeIDs) {
            my $repGenome2 = $repDB->rep_object($genome2);
            my $distance = $repGenome->distance($repGenome2);
            if ($distance < $min) { $min = $distance; }
            if ($distance > $max) { $max = $distance; }
            $sum += $distance;
            $count++;
            print $oh join("\t", $genome, $genome2, $distance) . "\n";
        }
        # Get the next genome.
        $genome = shift @genomeIDs;
        $num++;
    }
    close $oh;
    my $mean = $sum / $count;
    print STDERR "Distances range from $min to $max with a mean of $mean.\n";
}
