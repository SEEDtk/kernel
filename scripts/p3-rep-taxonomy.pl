=head1 Display Taxonomic Statistics on Representative Genomes

    p3-rep-taxonomy.pl [options] repDir

This script produces a report of the represented sets in a representative-genome database based on genus and species.
In each represented set, the number of genomes in each genus, and each species within that genus, will be output.
This is a long process because the taxonomic information is all in the PATRIC database and must be fetched.

Status messages will be written to the standard error output.

=head2 Parameters

The single positional parameter is the name of the directory containing the database. The files in the directory are
described in L<RepGenomeDb>.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use Math::Round;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir');
# Get the input directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Invalid or missing input directory $repDir.";
}
# Connect to PATRIC.
my $p3 = P3DataAPI->new();
# Load the database.
print STDERR "Loading database from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Get the list of genomes.
print STDERR "Counting set sizes.\n";
my $genomeList = $repDB->rep_list();
# Create a hash mapping each genome to its number of representatives.
my %repHash;
# This is a two-level hash with the genus totals.
my %genHash;
# This is a three-level hash with the genus and species information.
my %subHash;
my $counter = 0;
# Loop through the representative genomes.
for my $genomeID (@$genomeList) {
    my $repGenome = $repDB->rep_object($genomeID);
    # Get the represented genomes.  Note we add the representative itself.
    my $repList = $repGenome->rep_list();
    my @repGenomes = ($genomeID, map { $_->[0] } @$repList);
    # Create the species hash for this representative.
    my %gSubHash;
    # Get the genus and species information.
    $counter++;
    print STDERR "Retrieving taxonomy information for ($counter) " . $repGenome->name . ".";
    my $taxData = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name', 'genus', 'species'], \@repGenomes);
    print STDERR "  " . scalar(@$taxData) . " results.\n";
    for my $taxDatum (@$taxData) {
        # Extract the genus and species.  We fake it from the name if the main fields are blank.
        my ($genome, $name, $genus, $species) = @$taxDatum;
        my ($g2, $s2) = split /\s+/, $name;
        $genus ||= $g2; $species ||= "$g2 $s2";
        # Total up this genome.
        $repHash{$genomeID}++;
        $genHash{$genomeID}{$genus}++;
        $gSubHash{$genus}{$species}++;
    }
    # Save the taxonomic-breakdown hash.
    $subHash{$genomeID} = \%gSubHash;
}
# Sort the hash by represented-genome count.
print STDERR "Sorting groups.\n";
my @genomeIDs = sort { $repHash{$b} <=> $repHash{$a} } @$genomeList;
# We will accumulate useful statistics in here.
my $outliers = 0;
my %total = (genera => 0, species => 0, members => 0);
my $count = 0;
my $nCount = 0;
# Write the report.
print STDERR "Writing report.\n";
print join("\t", qw(id genus species count protLen name)) ."\n";
for my $genome (@genomeIDs) {
    $count++;
    my $repGenome = $repDB->rep_object($genome);
    my $protLen = length $repGenome->prot();
    my $name = $repGenome->name();
    my $genSubHash = $genHash{$genome};
    my $specSubHash = $subHash{$genome};
    my $specs = 0;
    for my $genus (keys %$specSubHash) {
        $specs += scalar(keys %{$specSubHash->{$genus}});
    }
    my $genCount = scalar keys %$genSubHash;
    my $memCount = $repHash{$genome};
    if ($memCount > 1) {
        $total{genera} += $genCount;
        $total{species} += $specs;
        $total{members} += $memCount;
        $nCount++;
    } else {
        $outliers++;
    }
    P3Utils::print_cols([$genome, $genCount, $specs, $memCount, $protLen, $name]);
}
# Now output the averages.
print STDERR "$count groups, $outliers singletons, $nCount real groups.\n";
print STDERR "Mean sizes for real groups: " . join(", ", map { nearest(0.1, $total{$_} / $nCount) . " $_" } qw(genera species members)) . ".\n";
print STDERR "All done.\n";