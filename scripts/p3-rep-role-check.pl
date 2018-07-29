=head1 Distance Validation Script

    p3-rep-role-check.pl [options]

This script checks the distance between genomes based on two measures-- kmer similarity in the seed protein, and role profile similarity.
The genome list is taken from the standard input. Every pair of genomes will be compared using the two measures. The seed protein name is
configurable, so that different seed proteins can be used; however, if a genome has more than one seed protein instance, this script does
not do any intelligent choosing, so a reliably universal role should be chosen.

The seed protein similarity is the number of protein kmers in common between the two seed protein instances.

The role profile similarity is the number of roles in common between the two genomes.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The L<P3Utils/col_options> can be used to choose the columns containing the genome ID.

The following additional options are supported.

=over 4

=item seedProt

The descriptive name of the seed protein to use for the protein kmer comparison. The default is C<Phenylalanyl tRNA-synthetase alpha chain>.

=item K

The kmer size to use for protein comparisons. The default is C<8>.

=item verbose

Write status messages to STDERR.

=item roleFile

The location of the C<roles.in.subsystems> file containing the roles of interest. The default is the one in the SEEDtk global data directory.

=back

=cut

use strict;
use FIG_Config;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use RoleParse;
use Stats;
use GEO;
use EvalCon;

# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::data_options(), P3Utils::col_options(), P3Utils::ih_options(),
        ['seedProt|seedprot|seed=s', 'seed protein description', { default => 'Phenylalanyl-tRNA synthetase alpha chain' }],
        ['kmer|K=i', 'protein kmer size', { default => 8 }],
        ['verbose|debug|v', 'show status messages on STDERR'],
        ['roleFile|rolefile|roles|r=s', 'roles.in.subsystems file', { default => "$FIG_Config::global/roles.in.subsystems" }]
        );
my $stats = Stats->new();
# Get the debug flag.
my $debug = $opt->verbose;
# Get access to PATRIC.
print STDERR "Connecting to database.\n" if $debug;
my $p3 = P3DataAPI->new();
# Get the checksum for the seed protein.
my $seedName = $opt->seedprot;
my $seedCheck = RoleParse::Checksum($seedName);
# Get the Kmer size.
my $K = $opt->kmer;
# Create an empty representative-genome set.
my $repDB = RepGenomeDb->new(K => $K);
# Read the role definitions.
print STDERR "Loading role definitions.\n" if $debug;
my ($nMap, $cMap) = EvalCon::LoadRoleHashes($opt->rolefile, $stats);
# Compute the options for GEO creation.
my $logH = ($debug ? \*STDERR : undef);
my %geoOptions = (roleHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, detail => 0, logH => $logH);
# This hash will map each genome pair {ID1\tID2} to its role profile similarity.
my %roleDist;
# This hash will map each genome pair {ID1\tID2} to its seed protein similarity.
my %seedDist;
# This hash contains the GEO for each genome.
my %geos;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($inHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Write the output headers.
if (! $opt->nohead) {
    P3Utils::print_cols(['genome1', 'name1', 'genome2', 'name2', 'kmerDiff', 'roleDiff']);
}
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # Get the seed protein for all the incoming genomes.
    print STDERR "Searching for seed proteins for " . scalar(@$couplets) . " genomes.\n" if $debug;
    my $protResults = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', $seedName]], ['genome_id', 'product', 'aa_sequence'], [map { $_->[0] } @$couplets], 'genome_id');
    # Get the genome protein sequences.
    my %gProts;
    for my $protResult (@$protResults) {
        my ($genomeID, $roleName, $prot) = @$protResult;
        if (RoleParse::Checksum($roleName) eq $seedCheck) {
            $gProts{$genomeID} = $prot;
        }
    }
    # Now get GEOs for these genomes.
    my @genomes = keys %gProts;
    print STDERR "Reading role profiles for " . scalar(@genomes) . " genomes.\n" if $debug;
    my $gHash = GEO->CreateFromPatric(\@genomes, %geoOptions);
    # Compare each genome to all the ones before it.
    for my $genome (@genomes) {
        my $geo = $gHash->{$genome};
        my $name = $geo->name;
        print STDERR "Processing $genome: $name.\n" if $debug;
        $stats->Add(genomeProcessed => 1);
        # Get the similarity scores for this genome. Specifying 0 for the score gets us all of them.
        my $gProt = $gProts{$genome};
        my $repList = $repDB->list_reps($gProt, 0);
        for my $genome2 (keys %geos) {
            # Compute the pair key.
            my $gKey = "$genome\t$genome2";
            # Store the protein similarity.
            $seedDist{$gKey} = $repList->{$genome2} // 0;
            # Compute the role similarity.
            my $geo2 = $geos{$genome2};
            $roleDist{$gKey} = $geo->role_similarity($geo2);
            $stats->Add(pairProcessed => 1);
        }
        # Add this genome to the databases.
        $repDB->AddRep($genome, $name, $gProt);
        $geos{$genome} = $geo;
    }
}
print STDERR "Sorting results.\n" if $debug;
my @pairs = sort { $seedDist{$b} <=> $seedDist{$a} } keys %seedDist;
for my $pair (@pairs) {
    my ($g1, $g2) = split /\t/, $pair;
    P3Utils::print_cols([$g1, $geos{$g1}->name, $g2, $geos{$g2}->name, $seedDist{$pair}, $roleDist{$pair}]);
    $stats->Add(pairOut => 1);
}
print STDERR "All done.\n" . $stats->Show() if $debug;
