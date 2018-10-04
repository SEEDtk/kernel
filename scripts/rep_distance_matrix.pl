=head1 Compute Distance Matrix for Representative Genomes

    rep_distance_matrix.pl [options] repDir

This script will take as input a representative-genome directory and compute the similarities between the genomes
represented. The resulting matrix will be output as a tab-delimited file. Each row will contain a genome ID and
then list the IDs of the other genomes, from closest to furthest, until a minimum similarity is reached (default C<25>).
The output file name will be C<repLists.tbl>. This file is used by the method for finding genomes in a region surrounding
a target genome-- L<RepGenomeDb/FindRegion>.

=head2 Parameters

The single positional parameter is the name of the directory containing the representative-genome database.

The command-line options are as follows.

=over 4

=item minScore

The minimum allowable similarity score. Representatives with a score equal to or higher than this will be included in the output.
The default is C<25>. Note that a score higher than the score used to build the representative set will result in no results for
any genome.

=back

=cut

use strict;
use P3Utils;
use RepGenomeDb;
use Stats;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir',
        ['minScore|min|m=i', 'minimum acceptable similarity score', { default => 25 }],
        );
my $stats = Stats->new();
# Check the parameters.
my $minScore = $opt->minscore;
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No representative genome directory specified.";
} elsif (! -d $repDir) {
    die "$repDir is either missing or invalid.";
}
# Load the rep-genome database.
print "Loading representative genomes from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir, verbose => 1, unconnected => 1);
# Get the list of genomes in the database.
my $genomes = $repDB->rep_list();
print scalar(@$genomes) . " representative genomes to process.\n";
# Open the output file.
open(my $oh, ">$repDir/repLists.tbl") || die "Could not open repLists.tbl: $!";
# Loop through the representative genomes.
for my $genome (@$genomes) {
    # Get this genome's seed protein.
    $stats->Add(genomesIn => 1);
    my $prot = $repDB->rep_object($genome)->prot;
    my $neighborH = $repDB->list_reps($prot, $minScore);
    # Get the genomes found in order. Note we have to pull the genome itself, since it is always the best match.
    my @found = sort { $neighborH->{$b} <=> $neighborH->{$a} } grep { $_ ne $genome } keys %$neighborH;
    my $found = scalar @found;
    print "$found close reps found for $genome.\n";
    $stats->Add(neighborsFound => $found);
    P3Utils::print_cols([$genome, @found], oh => $oh);
}
print "All done.\n" . $stats->Show();
