=head1 Simance Validation Script

    p3-rep-role-check.pl [options] gtoDir

This script checks the Simance between genomes based on two measures-- kmer similarity in the seed protein, and role profile similarity.
The genome list is taken from the standard input. Every pair of genomes will be compared using the two measures. The seed protein name is
configurable, so that different seed proteins can be used; however, if a genome has more than one seed protein instance, this script does
not do any intelligent choosing; therefore, a reliably universal role should be chosen.

The seed protein similarity is the number of protein kmers in common between the two seed protein instances.

The role profile similarity is the number of roles in common between the two genomes divided by the number of roles found in either genome.

=head2 Parameters

The positional parameter is the name of a directory containing L<GenomeTypeObject> files (all with the filename suffix C<.gto>). These will
be the genomes being compared.

The command-line options in L<EvalCon::role_options> are used to select the role definition files. The following additional options are supported.

=over 4

=item seedProt

The ID of the seed protein to use for the protein kmer comparison. The default is C<PhenTrnaSyntAlph>.

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
use P3Utils;
use RepGenomeDb;
use RoleParse;
use Stats;
use GEO;
use EvalCon;

# Get the command-line options.
my $opt = P3Utils::script_opts('inDir', P3Utils::data_options(), P3Utils::col_options(), P3Utils::ih_options(),
        ['seedProt|seedprot|seed=s', 'seed protein description', { default => 'PhenTrnaSyntAlph' }],
        ['kmer|K=i', 'protein kmer size', { default => 8 }],
        ['verbose|debug|v', 'show status messages on STDERR'],
        EvalCon::role_options(),
        );
my $stats = Stats->new();
# Get the debug flag.
my $debug = $opt->verbose;
my $logH = ($debug ? \*STDERR : undef);
# Check the parameters.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Missing or invalid input directory $inDir.";
}
# Read the input directory.
print STDERR "Reading $inDir.\n" if $debug;
opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
my @gtos = grep { $_ =~ /\S\.gto$/ } readdir $dh;
closedir $dh;
print STDERR scalar(@gtos) . " input files found in $inDir.\n" if $debug;
# Get the checksum for the seed protein.
my $seedName = $opt->seedprot;
# Get the Kmer size.
my $K = $opt->kmer;
# Create an empty representative-genome set.
my $repDB = RepGenomeDb->new(K => $K);
# Read the role definitions.
my $evalCon = EvalCon->new_for_script($opt, $logH);
my ($nMap, $cMap) = $evalCon->roleHashes;
my %rolesToUseH = map { $_ => 1 } @{$evalCon->rolesToUse};
# Compute the options for GEO creation.
my %geoOptions = (roleHashes => [$nMap, $cMap], stats => $stats, detail => 2, logH => $logH, rolesToUse => \%rolesToUseH);
# This hash will map each genome pair {ID1\tID2} to its role profile similarity.
my %roleSim;
# This hash will map each genome pair {ID1\tID2} to its seed protein similarity.
my %seedSim;
# This hash contains the GEO for each genome.
my %geos;
# Write the output headers.
P3Utils::print_cols(['genome1', 'name1', 'genome2', 'name2', 'kmerDiff', 'roleDiff']);
# Loop through the input directory.
for my $gto (@gtos) {
    # Create the GTO.
    my $gHash = GEO->CreateFromGtoFiles(["$inDir/$gto"], %geoOptions);
    my ($genome) = keys %$gHash;
    die "Error loading from $gto." if ! $genome;
    my $geo = $gHash->{$genome};
    # Get the seed protein.
    my $fidList = $geo->roleFids($seedName);
    if (! @$fidList) {
        die "No seed protein found in $genome.";
    }
    my $gProt = $geo->protein($fidList->[0]);
    if (! $gProt) {
        die "Invalid seed protein in $genome.";
    }
    my $name = $geo->name;
    print STDERR "Processing $genome: $name.\n" if $debug;
    $stats->Add(genomeProcessed => 1);
    # Get the similarity scores for this genome. Specifying 0 for the score gets us all of them.
    my $repList = $repDB->list_reps($gProt, 0);
    for my $genome2 (keys %geos) {
        # Compute the pair key.
        my $gKey = "$genome\t$genome2";
        # Store the protein similarity.
        $seedSim{$gKey} = $repList->{$genome2} // 0;
        # Compute the role similarity.
        my $geo2 = $geos{$genome2};
        $roleSim{$gKey} = $geo->role_similarity($geo2);
        $stats->Add(pairProcessed => 1);
    }
    # Add this genome to the databases.
    $repDB->AddRep($genome, $name, $gProt);
    $geos{$genome} = $geo;
}
print STDERR "Sorting results.\n" if $debug;
my @pairs = sort { $seedSim{$b} <=> $seedSim{$a} } keys %seedSim;
for my $pair (@pairs) {
    my ($g1, $g2) = split /\t/, $pair;
    P3Utils::print_cols([$g1, $geos{$g1}->name, $g2, $geos{$g2}->name, $seedSim{$pair}, $roleSim{$pair}]);
    $stats->Add(pairOut => 1);
}
print STDERR "All done.\n" . $stats->Show() if $debug;
