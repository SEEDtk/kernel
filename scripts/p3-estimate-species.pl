=head1 Estimate Species of Contigs in a FASTA File

    p3-estimate-species.pl [options] fastaFile

This script searches a FASTA file for a key protein and uses it to compute the species of a genome described by a FASTA
file.  The default is to use the universal protein Phenylalanyl tRNA-synthetase alpha chain. There is almost always exactly
one of these per genome.  If we do not find a good candidate for the protein or find more than one good candidate, the
script will fail.

The standard output will a single tab-delimited line containing the taxonomic ID of the species, the name of the species,
and then the full taxonomy, names delimited by semi-colons, from the domain down to the species group.

=head2 Parameters

The single positional parameter is the name of a FASTA file containing the contigs of the proposed genome.

The command-line options are as follows.

=over 4

=item seedProtFasta

A FASTA file containing examples of the universal role to use for seeding the bin assignment.  The default is
C<seedprot.fa> in the global data directory.

=item seedfasta

The name of the BLAST database for the seed protein in the various PATRIC genomes. The default is
C<PhenTrnaSyntAlph.fa> in the global data directory.

=item maxE

The maximum acceptable E-value. The default is C<1e-20>.

=item refMaxE

The maximum acceptable E-value for blasting to determine the best reference genome for a seed contig. Each seed
contig eventually forms a bin. The default is C<1e-10>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=item verbose

Write status messages to STDERR.

=item rank

The desired accuracy level. This should be a taxonomic rank. The default is C<species>. Other values are
(for example) C<genus> and C<strain>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Bin::Blast;
use Loader;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('fastaFile',
                ['seedProtFasta=s', 'name of a FASTA file containing examples of the seed protein to locate',
                                    { default => "$FIG_Config::global/seedprot.fa" }],
                ['seedfasta|F=s',  'BLAST database (or FASTA file) of seed protein in all genomes', { default => "$FIG_Config::global/PhenTrnaSyntAlph.fa"}],
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-20 }],
                ['refMaxE=f',      'maximum acceptable e-value for reference genome blast hits', { default => 1e-10 }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
                ['verbose|debug|v', 'write status messages to STDERR'],
                ['rank=s',         'rank level of desired output', { default => 'species' }],
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Check the parameters.
my ($fastaFile) = @ARGV;
if (! $fastaFile) {
    die "No input file specified.";
} elsif (! -s $fastaFile) {
    die "Input file $fastaFile missing or empty.";
}
# Compute the rank.
my $rank = $opt->rank;
# Get the debug flag.
my $debug = $opt->verbose;
# Compute the blast parameters.
my $maxE = $opt->maxe;
my $rMaxE = $opt->refmaxe // $maxE;
# Get the I/O utility object.
my $loader = Loader->new();
# Create the blaster.
my $blaster = Bin::Blast->new(undef, $FIG_Config::temp, $fastaFile,
        maxE => $opt->maxe, minlen => $opt->minlen, gap => $opt->gap, silent => 1);
# Find the seed protein.
my $protFile = $opt->seedprotfasta;
print STDERR "Searching for seed protein in $fastaFile using $protFile.\n" if $debug;
my $matches = $blaster->FindProtein($protFile);
# Only proceed if we found a single match.
my ($contig, @others) = keys %$matches;
if (! $contig) {
    print STDERR "No seed proteins found.\n" if $debug;
} elsif (scalar @others) {
    print STDERR "Too many seed proteins found.\n" if $debug;
} else {
    # Get the DNA for the protein.
    my $seqHash = $loader->GetDNA($matches, $fastaFile);
    # Get the best match for the DNA.
    my $seedFastaFile = $opt->seedfasta;
    print STDERR "Searching for best DNA match to seed protein in $contig using $seedFastaFile.\n" if $debug;
    my $contigHash = $blaster->MatchProteins($seqHash, undef, 1, $opt->refmaxe, db => $seedFastaFile, type => 'dna');
    # MatchProteins is a general-purpose thing that returns a list of genomes for each contig. The "1" in the parameter
    # list insures that it returns at most one.
    my $genomeData = $contigHash->{$contig}[0];
    if (! $genomeData) {
        print STDERR "No matching genome found for $contig.\n" if $debug;
    } else {
        my ($genome, $score, $name) = @$genomeData;
        print STDERR "Retrieving taxonomic lineage from $genome.\n" if $debug;
        my $genomeRecords = P3Utils::get_data_keyed($p3, genome => [], ['taxon_lineage_ids'], [$genome]);
        if (! @$genomeRecords) {
            print STDERR "Genome $genome not found in PATRIC.\n" if $debug;
        } else {
            my $taxonList = $genomeRecords->[0][0];
            if (! $taxonList) {
                print STDERR "Lineage not available for $genome.\n" if $debug;
            } else {
                print STDERR "Looking for $rank in taxonomy list.\n" if $debug;
                my $taxonRecords = P3Utils::get_data_keyed($p3, taxonomy => [['eq', 'taxon_rank', $rank]],
                        ['taxon_id', 'taxon_name', 'lineage_names'], $taxonList);
                if (! @$taxonRecords) {
                    print STDERR "No $rank found in lineage for $genome.\n" if $debug;
                } else {
                    my $taxonRecord = $taxonRecords->[0];
                    my ($id, $name, $lineage) = @$taxonRecord;
                    if (! $lineage) {
                        print STDERR "No lineage found for species $id.\n" if $debug;
                    } else {
                        # Strip off any super-domains.
                        while ($lineage->[0] =~ /^(?:root|cellular)/i) {
                            shift @$lineage;
                        }
                        # Write the result.
                        print join("\t", $id, $name, join('; ', @$lineage)) . "\n";
                    }
                }
            }
        }
    }
}