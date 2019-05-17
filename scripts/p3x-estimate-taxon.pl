=head1 Estimate Taxonomic Grouping of Contigs in a FASTA File

    p3-estimate-taxon.pl [options] fastaFile1 fastaFile2 ... fastaFileN

This script searches a FASTA file for a key protein and uses it to compute the species of a genome described by a FASTA
file.  The default is to use the universal protein Phenylalanyl tRNA-synthetase alpha chain. There is almost always exactly
one of these per genome.  If we do not find a good candidate for the protein or find more than one good candidate, the
script will fail.

The standard output will a single tab-delimited line containing the input file name, the taxonomic ID, the name,
and then the full taxonomy, names delimited by semi-colons, from the domain down to the species group.

=head2 Parameters

The positional parameters are the names of FASTA files for the proposed genomes.

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

=item input

If specified, then should be the name of a file containing the names of FASTA files in the first column.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use TaxCheck;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('fastaFile1 fastaFile2 ... fastaFileN',
                ['seedProtFasta=s', 'name of a FASTA file containing examples of the seed protein to locate',
                                    { default => "$FIG_Config::p3data/seedprot.fa" }],
                ['seedfasta|F=s', 'BLAST database (or FASTA file) of seed protein in all genomes', { default => "$FIG_Config::p3data/PhenTrnaSyntAlph.fa"}],
                ['maxE|e=f', 'maximum acceptable e-value for blast hits', { default => 1e-20 }],
                ['refMaxE=f', 'maximum acceptable e-value for reference genome blast hits', { default => 1e-10 }],
                ['gap|g=i', 'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['minlen|l=f', 'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
                ['verbose|debug|v', 'write status messages to STDERR'],
                ['rank=s', 'rank level of desired output', { default => 'genus' }],
                ['input=s', 'optional input file']
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Check the parameters.
my (@fastaFiles) = @ARGV;
if (! @fastaFiles && ! $opt->input) {
    die "No input files specified.";
}
# Compute the rank.
my $rank = $opt->rank;
# Get the debug flag.
my $debug = $opt->verbose;
# Compute the blast parameters.
my $maxE = $opt->maxe;
my $rMaxE = $opt->refmaxe;
my $gap = $opt->gap;
my $minlen = $opt->minlen;
# Get the protein files.
my $seedFastaFile = $opt->seedfasta;
my $protFile = $opt->seedprotfasta;
# Create the taxonomy checker.
my $taxChecker = TaxCheck->new($p3, debug => $debug, protFile => $protFile, seedFastaFile => $seedFastaFile,
        maxE => $maxE, refMaxE => $rMaxE, gap => $gap, minlen => $minlen);
# Check for additional input.
if ($opt->input) {
    open(my $ih, '<', $opt->input) || die "Could not open input file: $!";
    my $fileColumn = P3Utils::get_col($ih, 0);
    push @fastaFiles, @$fileColumn;
    print STDERR scalar(@$fileColumn) . " file names read from input.\n" if $debug;
}
# Loop through the FASTA files, computing the taxonomic groups.
for my $fastaFile (@fastaFiles) {
    if (! -s $fastaFile) {
        print STDERR "$fastaFile not found.\n" if $debug;
    } else {
        my ($id, $name, $lineage) = $taxChecker->Compute($fastaFile, $rank);
        if ($id) {
            my $lineageString = join('; ', @$lineage);
            P3Utils::print_cols([$id, $name, $lineageString]);
        }
    }
}