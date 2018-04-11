=head1 Filter Genomes For Good Seed Proteins

    p3-good-seed.pl [options]

This script takes as input a set of genome IDs and filters them for the presence of a good seed protein. The seed protein
is Phenylalanyl tRNA-synthetase alpha chain, a protein that occurs singly in virtually all genomes. We look for all
occurrences. If there is more than one, or its length is inappropriate to the domain, we reject the genome. Otherwise,
the input line is copied to the output unchanged.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The key column (containing the genome IDs) can be specified using L<P3Utils/col_options>.
=cut

use strict;
use P3DataAPI;
use P3Utils;
use RoleParse;
use SeedUtils;


# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::col_options(), P3Utils::ih_options(),
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Form the full header set and write it out.
if (! $opt->nohead) {
    P3Utils::print_cols($outHeaders);
}
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # The first task is to get the taxonomic information so we can compute the domains.
    my @genomeIDs = map { $_->[0] } @$couplets;
    my $resultList = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'taxon_lineage_names'], \@genomeIDs, 'genome_id');
    # Build a hash of genome IDs to domains.
    my %gDomain;
    for my $result (@$resultList) {
        my ($genomeID, $taxons) = @$result;
        my $domain = $taxons->[1];
        if ( !grep { $_ eq $domain } qw(Bacteria Archaea Eukaryota) ) {
            $domain = $taxons->[0];
        }
        $gDomain{$genomeID} = $domain;
    }
    # Now we analyze the seed proteins. Note that we have to ask for the amino acid sequence. The length field does not
    # always work.
    $resultList = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']],
            ['genome_id', 'product', 'aa_sequence'], \@genomeIDs, 'genome_id');
    # Loop through the results, computing good and bad. This hash will map each genome ID to the number of seed proteins.
    # A bad seed protein counts as 2. We consider a genome good if its ID is mapped to 1.
    my %good;
    for my $result (@$resultList) {
        my ($genomeID, $function, $protein) = @$result;
        my $domain = $gDomain{$genomeID};
        # Genome is only good if we know the domain!
        if ($domain) {
            # Now we check the function.
            my $checksum = RoleParse::Checksum($function);
            if ($checksum eq 'WCzieTC/aZ6262l19bwqgw') {
                my ($min, $max) = (209, 405);
                if ($domain eq 'Archaea') {
                    ($min, $max) = (293, 652);
                }
                my $aaLen = length $protein;
                if ($aaLen >= $min && $aaLen <= $max) {
                    $good{$genomeID} += 1;
                } else {
                    $good{$genomeID} += 2;
                }
            }
        }
    }
    # Now we run through the couplets and output any couplet for which $good{$genomeID} == 1.
    for my $couplet (@$couplets) {
        my ($genomeID, $line) = @$couplet;
        if ($good{$genomeID} && $good{$genomeID} == 1) {
            P3Utils::print_cols($line);
        }
    }
}