=head1 Find the Neighborhood of a Genome

    find_neighborhood.pl [options] repDir genomeID

This script will find the genomes believed to be in the neighborhood of another genome using the data in a representative-genomes
directory. The output will contain the neighboring genomes and their similarity scores compared to the original in the form
of a two-column tab-delimited file.

=head2 Parameters

The positional parameters are the name of the representative-genome directory, which must contain the files described in
L<RepGenomeDb/new_from_dir> and the ID of the genome whose neighborhood is desired. If the genome has a missing seed protein
no neighborhood will be returned. If the genome is in a sparse part of the phylogenetic tree the return set will be small.

The command-line options are as follows.

=over 4

=item size

The desired size of the neighborhood. The default is C<100>.

=item minScore

The minimum similarity score for a representative genome for its represented set to be considered close. The default is C<25>.

=item sliceSize

The number of genomes to take from each representative set. The default is C<35>.

=item verbose

If specified, status messages will be written to STDERR.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RoleParse;
use RepGenomeDb;
use RepGenome;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir genomeID', P3Utils::data_options(), P3Utils::col_options(), P3Utils::ih_options(),
        ['size', 'neighborhood size', { default => 100 }],
        ['minScore', 'minimum similarity score for representatives', { default => 25 }],
        ['sliceSize', 'maximum number of genomes per represented set', { default => 35 }],
        ['verbose|debug|v', 'write status messages to STDERR'],
        );
my $debug = $opt->verbose;
# Get the parameters.
my ($repDir, $genomeID) = @ARGV;
if (! $repDir) {
    die "No representative-genome directory specified.";
} elsif (! -d $repDir) {
    die "Representative-genome directory $repDir is missing or invalid.";
} elsif (! $genomeID) {
    die "No input genome ID specified.";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Write the output file header.
P3Utils::print_cols(['genomeID', 'score']);
# Get the seed protein of the input genome.
print STDERR "Locating $genomeID seed protein.\n" if $debug;
my $seedHash = FindSeeds([$genomeID]);
my $gProt = $seedHash->{$genomeID};
if (! $gProt) {
    die "No seed protein found in $genomeID.";
}
# Load the representative-genomes database.
print STDERR "Loading representative genomes from $repDir.\n" if $debug;
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Find the neighbors.
print STDERR "Searching for neighbors.\n" if $debug;
my $neighbors = $repDB->FindRegion($gProt, size => $opt->size, sliceSize => $opt->slicesize, minScore => $opt->minscore);
print STDERR scalar(@$neighbors) . " neighboring genomes found.\n" if $debug;
# Get their seed proteins.
$seedHash = FindSeeds($neighbors);
# Create a rep-genome object for the main protein.
my $repObject = RepGenome->new($genomeID, prot => $gProt, K => $repDB->K);
# Get the similarities for the other genomes and print them.
for my $other (@$neighbors) {
    my $prot = $seedHash->{$other};
    die "No protein found for $other." if ! $prot;
    my $score = $repObject->check_genome($prot);
    P3Utils::print_cols([$other, $score]);
}

## Get the seed proteins of the listed genomes. Returns a hash ref.
sub FindSeeds {
    my ($genomes) = @_;
    # Get all the features that look promising.
    my $results = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']],
            ['genome_id', 'product', 'aa_sequence'], $genomes, 'genome_id');
    # The ones we keep go in here.
    my %retVal;
    # Loop through the results.
    for my $genomeDatum (@$results) {
        my ($id, $function, $prot) = @$genomeDatum;
        my $checksum = RoleParse::Checksum($function);
        if ($checksum eq 'WCzieTC/aZ6262l19bwqgw') {
            $retVal{$id} = $prot;
        }
    }
    # Return the proteins found.
    return \%retVal;
}
