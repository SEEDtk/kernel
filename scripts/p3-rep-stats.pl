=head1 Display Statistics on Representative Genomes

    p3-rep-stats.pl [options] repDir

This script produces a simple report on the number of genomes represented by each genome in a representative-genome database.

=head2 Parameters

The single positional parameter is the name of the directory containing the database. The files in the directory are
described in L<RepGenomeDb>.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir');
# Get the input directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Invalid or missing input directory $repDir.";
}
# Load the database.
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Get the list of genomes.
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
print join("\t", qw(id count protLen name)) ."\n";
for my $genome (@genomeIDs) {
    my $repGenome = $repDB->rep_object($genome);
    my $protLen = length $repGenome->prot();
    my $name = $repGenome->name();
    print join("\t", $genome, $repHash{$genome}, $protLen, $name) . "\n";
}

