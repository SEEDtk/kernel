=head1 Build Reference-Genome Table from Sorted Evaluations

    p3-eval-ref-table.pl [options] outDir

This script builds the reference-genome table from BV-BRC for each genus and species for which a good genome can be found.
The reference-genome table is a tab-delimited file, each record containing (0) a taxonomic grouping ID and (1) the best-quality
genome found in that grouping. The groupings are all species or genus level. The table is produced under the name C<ref.genomes.tbl>
in the specified output directory.

The input file must have been produced by L<p3x-eval-sort.pl>. This script will be faster if there is a C<taxon_lineage_ids>.

=head2 Parameters

The positional parameter is the name of the output directory to contain the reference genomes.

The standard input can be overridden using the options in L<P3Utils/ih_options>. The database can be specified using the options in
L<Shrub/script_options>.

The following additional command-line options are supported.

=over 4

=item batchSize

Maximum number of lines to read in a batch. The default is C<1000>.

=item merge

If specified, the name of a file containing the merge data from the latest NCBI taxonomy dump.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Shrub;
use Stats;
use Data::Dumper;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir', Shrub::script_options(), P3Utils::ih_options(),
                ['batchSize|b=i', 'input batch size', { default => 1000 }],
                ['merge=s', 'NCBI merged.dmp file name']
        );
my $stats = Stats->new();
# Check the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Output directory $outDir missing or invalid.";
}
print "Connecting to databases.\n";
# Get access to BV-BRC.
my $p3 = P3DataAPI->new();
# Get access to the Shrub.
my $shrub = Shrub->new_for_script($opt);
# Create hashes for each genus and species.
my %genus = map { $_ => undef } $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(type) = ?', ['genus'], 'id');
print scalar(keys %genus) . " genus groupings found in Shrub.\n";
my %species = map { $_ => undef } $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(type) = ?', ['species'], 'id');
print scalar(keys %species) . " species groupings found in Shrub.\n";
# Get the merge file. We need to map each old taxon to its new version.
if ($opt->merge) {
    open(my $mh, '<', $opt->merge) || die "Could not open merged.dmp: $!";
    while (! eof $mh) {
        my $line = <$mh>;
        if ($line =~ /(\d+)\t\|\t(\d+)/) {
            my ($old, $new) = ($1, $2);
            if ($species{$new}) {
                $species{$old} = undef;
            } elsif ($genus{$new}) {
                $genus{$new} = undef;
            }
            $stats->Add(taxonMergeRead => 1);
        }
    }
}
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $cols) = P3Utils::find_headers($ih, qualityFile => 'genome_id', 'Good Genome');
my ($keyCol, $goodFlagCol) = @$cols;
# Check for a taxon lineage column.
my $taxCol = P3Utils::find_column('taxon_lineage_ids', $outHeaders, 'optional');
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    my @goodGenomes = map { $_->[0] } grep { $_->[1][$goodFlagCol] } @$couplets;
    if (scalar(@goodGenomes)) {
        print scalar(@goodGenomes) . " good genomes found in batch.\n";
        # Get the lineage for each group.
        my $taxList;
        if (defined $taxCol) {
            $taxList = [];
            for my $couplet (@$couplets) {
                push @$taxList, [$couplet->[0], [split /::/, $couplet->[1][$taxCol]]];
            }
        } else {
            $taxList = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'taxon_lineage_ids'], \@goodGenomes);
        }
        # Form it into a hash. Note that sometimes the lineage comes back as an empty string instead of a list, so
        # we do the or-thing.
        my %gTaxes = map { $_->[0] => ($_->[1] || []) } @$taxList;
        # Now loop through the lineage for each genome. If the genus and species do not already have reference genomes, save them.
        for my $genome (@goodGenomes) {
            $stats->Add(goodGenomeChecked => 1);
            my $taxList = $gTaxes{$genome};
            if (! $taxList) {
                $stats->Add(lineageNotFound => 1);
            } else {
                my $used;
                for my $taxon (@$taxList) {
                    $stats->Add(taxonChecked => 1);
                    if (exists $species{$taxon} && ! $species{$taxon}) {
                        $species{$taxon} = $genome;
                        $stats->Add(speciesReference => 1);
                        $used = 1;
                    } elsif (exists $genus{$taxon} && ! $genus{$taxon}) {
                        $genus{$taxon} = $genome;
                        $stats->Add(genusReference => 1);
                        $used = 1;
                    }
                }
                if ($used) {
                    print "$genome chosen as a reference.\n";
                    $stats->Add(refGenomeFound => 1);
                }
            }
        }
    }
}
print "Writing output.\n";
open(my $oh, ">$outDir/ref.genomes.tbl") || die "Could not open ref.genomes.tbl: $!";
for my $hash (\%genus, \%species) {
    for my $taxon (sort keys %$hash) {
        my $genome = $hash->{$taxon};
        if ($genome) {
            print $oh "$taxon\t$genome\n";
        }
    }
}
print "All done\n" . $stats->Show();