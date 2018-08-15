=head1 Create Neighbor Sets

    create_neighbor_sets.pl [options] repDir outDir

This script takes as input a list of genome IDs and a representative-genomes directory. Each genome should be from a different
genus. For each genome, the script produces a directory with the genus name and containing a C<rep.genomes> file with the IDs
of neighboring genomes. Each directory can be fed into the What's-Changed pipeline to create a tree structure.

=head2 Parameters

The positional parameters are the name of the representative-genomes directory and the name of the output directory. The output
will be in subdirectories created under the output directory.

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain genome IDs in the key column.
The key column can be specified using the options in L<P3Utils/col_options>.

The following additional options are supported.

=over 4

=item clear

If specified, the output directory will be erased before beginning.

=item size

The desired size of the neighborhood for each input genome. The default is C<100>.

=item slice

The size of a slice. The neighborhood is divided into slices at progressively further distances. This option specifies the maximum slice
size.  The default is C<35>.

=item min

The minimum acceptable score for a represented set to be considered a neighbor. The default is C<25>.

=item seed

The name of the seed protein. This should be a protein description string with no EC or TC number. The default is
C<Phenylalanyl-tRNA synthetase alpha chain>. The protein must match the one used by the representative-genomes
directory.

=item minSize

The minimum size for a neighborhood to be considered significant.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use RoleParse;
use Stats;

$| = 1;
my $stats = Stats->new();
# Get the command-line options.
my $opt = P3Utils::script_opts('repDir outDir', P3Utils::col_options(), P3Utils::ih_options(),
        ['min|m=i', 'minimum similarity score', { default => 25 }],
        ['slice|S=i', 'maximum slice size', { default => 35 }],
        ['size|s=i', 'optimal neighborhood size', { default => 100 }],
        ['clear', 'clear output directory before beginning'],
        ['seed=s', 'name of the seed protein to use', { default => 'Phenylalanyl-tRNA synthetase alpha chain' }],
        ['minSize=i', 'minimum neighborhood size', { default => 50 }]
        );
# Get the input parameters.
my ($repDir, $outDir) = @ARGV;
if (! $repDir) {
    die "No representative-genome directory specified."
} elsif (! -s "$repDir/6.1.1.20.fasta") {
    die "$repDir does not appear to be a representative genomes directory.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} elsif ($opt->clear) {
    print "Erasing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase $outDir: $!";
}
# Load the RepGenome database.
print "Loading from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir, verbose => 1);
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Compute the name and checksum of the seed protein.
my $protName = $opt->seed;
my $protCheck = RoleParse::Checksum($protName);
print "Seed protein is $protName.\n";
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    $stats->Add(batchIn => 1);
    # Get the name and genus for each genome.
    my @genomes = map { $_->[0] } @$couplets;
    my $count = scalar(@genomes);
    $stats->Add(genomeIn => $count);
    print "Extracting genome data for $count genomes.\n";
    my $gResults = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name', 'genus'], \@genomes);
    print "Extracting protein data for $count genomes.\n";
    my $pResults = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', $protName]], ['genome_id', 'product', 'aa_sequence'], \@genomes, 'genome_id');
    # Form the genome data into a hash.
    my %gResults = map { shift(@$_) => $_ } @$gResults;
    # Loop through the proteins. We will append the AA sequence to the end of each genome's data list.
    print "Analyzing proteins.\n";
    for my $pResult (@$pResults) {
        $stats->Add(seedIn => 1);
        my ($genome, $product, $prot) = @$pResult;
        # Validate the checksum.
        my $check = RoleParse::Checksum($product);
        if ($check ne $protCheck) {
            $stats->Add(seedFunny => 1);
        } else {
            my $gDatum = $gResults{$genome};
            if (! $gDatum) {
                print "Skipping $genome: no genome data found.\n";
                $stats->Add(genomeNoData => 1);
            } elsif (! $gDatum->[1]) {
                print "Skipping $genome: no genus found.\n";
                $stats->Add(genomeBadGenus => 1);
            } else {
                push @$gDatum, $prot;
                $stats->Add(seedFound => 1);
            }
        }
    }
    # Now %gResult has [name, genus, protein] for each valid genome.
    for my $genome (sort keys %gResults) {
        my ($name, $genus, $prot) = @{$gResults{$genome}};
        if (! $prot) {
            print "Skipping $genome: no seed protein found.\n";
            $stats->Add(genomeNoSeed => 1);
        } else {
            $genus =~ s/\s+/_/g;
            my $outputDir = "$outDir/$genus";
            if (-d $outputDir) {
                print "Skipping $genome: output directory for $genus already exists.\n";
                $stats->Add(genomeRedundant => 1);
            } else {
                # Try to create a neighborhood.
                print "Processing $genome.\n";
                my $genomeList = $repDB->FindRegion($prot, size => $opt->size, sliceSize => $opt->slice, minScore => $opt->min);
                my $retCount = scalar @$genomeList;
                # Verify that it's worth writing to output.
                if ($retCount < $opt->minsize) {
                    print "Only $retCount neighbors found for $genome: no output.\n";
                    $stats->Add(neighborhoodSmall => 1);
                } else {
                    File::Copy::Recursive::pathmk($outputDir) || die "Could not create $outputDir: $!";
                    open(my $oh, ">$outputDir/rep.genomes") || die "Could not open output file in $outputDir: $!";
                    print $oh map { "$_\n" } @$genomeList;
                    $stats->Add(neighborsFound => scalar @$genomeList);
                    $stats->Add(neighborhoodsOut => 1);
                }
            }
        }
    }
}
print "All done.\n" . $stats->Show();