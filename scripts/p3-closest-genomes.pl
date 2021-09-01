=head1 Find the Closest Genomes Using a Representative-Genomes Directory

    p3-closest-genomes.pl [options] repDir genomeID

Use a representative-genomes directory to find the genomes closest to a given genome.  The standard output will be tab-delimited
and contain the ID, name, and similarity score for each genome found.

=head2 Parameters

The positional parameters are the name of a representative genome directory (see L<RepGenomeDb>) and the ID of a genome
whose neighbors are desired.

The following command-line parameters are supported.

=over 4

=item verbose

If specified, progress messages will be written to the standard error output.

=item role

The role for the seed protein.  The default is C<Phenylalanyl-tRNA synthetase alpha chain>.  The role chosen should match the
one used to build the representative-genome database.

=item N

The number of genomes to return.  The default is C<10>.

=item minScore

The minimum score for a neighboring representative set.  The default is C<25>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use P3RepGenomes;
use RepGenomeDb;


# Get the command-line options.
my $opt = P3Utils::script_opts('repDir genome',
        ['verbose|debug|v', 'write progress messages to the standard error output'],
        ['role|R=s', 'seed protein role', { default => 'Phenylalanyl-tRNA synthetase alpha chain' }],
        ['num|N=i', 'number of genomes to return', { default => 10 }],
        ['minScore|minscore|min=i', 'minimum score for neighboring sets', { default => 25 }]
        );
# Get the options.
my $debug = $opt->verbose;
my $role = $opt->role;
my $N = $opt->num;
my $minScore = $opt->minscore;
# Get access to BV-BRC.
print STDERR "Connecting to BV-BRC.\n" if $debug;
my $p3 = P3DataAPI->new();
# Get the positional parameters.
my ($repDir, $genome) = @ARGV;
if (! $repDir) {
    die "No representative-genomes directory specified.";
} elsif (! -d "$repDir") {
    die "Directory $repDir not found or invalid.";
} elsif (! $genome) {
    die "No genome ID specified.";
}
# Read the rep-genome directory.
print STDERR "Loading $repDir.\n" if $debug;
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Create the master object.
my $p3RepDB = P3RepGenomes->new($p3, $repDB, role => $role);
# We need the name and the protein sequence for the target genome.
my $gHash = $p3RepDB->GetProts([$genome]);
my $gData = $gHash->{$genome};
if (! $gData) {
    die "$genome not found in BV-BRC.";
} else {
    my ($gName, $seedProt) = @$gData;
    # Ask for the closest genomes.
    print STDERR "Looking for neighbors of $genome: $gName.\n";
    $gHash = $p3RepDB->FindClosest($seedProt, minScore => $minScore, N => $N);
    # Write the output headers.
    P3Utils::print_cols([qw(genome name score)]);
    # Compute the high score.
    my $highScore = length($seedProt) - $repDB->K;
    # Add the base genome to the list and sort it.
    print STDERR "Sorting results.\n" if $debug;
    $gHash->{$genome} = [$highScore, $gName];
    my @neighbors = sort { $gHash->{$b}[0] <=> $gHash->{$a}[0] } keys %$gHash;
    for my $nGenome (@neighbors) {
        $gData = $gHash->{$nGenome};
        P3Utils::print_cols([$nGenome, @$gData]);
    }
}
