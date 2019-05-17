use strict;
use P3DataAPI;
use P3Utils;
use RoleParse;
use SeedUtils;
use File::Copy::Recursive;
use Stats;


# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::col_options(), P3Utils::ih_options(),
        );
my $stats = Stats->new();
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Open the output file.
open(my $oh, ">CheckG/newMatrix.tbl") || die "Could not open output: $!";
my $headers = <$ih>;
print $oh $headers;
if (! -d 'TempSciKit') {
    File::Copy::Recursive::pathmk('TempSciKit') || die "Could not create temp dir: $!";
}
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, 0, $opt);
    # The first task is to get the taxonomic information so we can compute the domains.
    my @genomeIDs = map { $_->[0] } @$couplets;
    my $resultList = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'taxon_lineage_names'], \@genomeIDs, 'genome_id');
    # We consider each genome as good unless it is archaea. In that case, we have to do a consistency check.
    my %good;
    for my $result (@$resultList) {
        my ($genomeID, $taxons) = @$result;
        if (! $taxons) {
            print "No taxon lineage found for $genomeID.\n";
            $stats->Add(badTaxon => 1);
        } else {
            my $domain = $taxons->[1];
            if ( !grep { $_ eq $domain } qw(Bacteria Archaea Eukaryota) ) {
                $domain = $taxons->[0];
            }
            if ($domain ne 'Archaea') {
                $good{$genomeID} = 1;
                $stats->Add(nonArchaea => 1);
            } else {
                File::Copy::Recursive::pathempty('TempSciKit');
                print "Checking $genomeID.\n";
                my $gto = $p3->gto_of($genomeID);
                $gto->destroy_to_file('TempSciKit/bin.gto');
                # Compute the quality.
                my $cmd = "gto_consistency TempSciKit/bin.gto TempSciKit/results $FIG_Config::p3data/FunctionPredictors $FIG_Config::p3data/roles.in.subsystems $FIG_Config::p3data/roles.to.use";
                SeedUtils::run($cmd);
                # Read in the results.
                open(my $qh, "<TempSciKit/results/evaluate.log") || die "Could not open $genomeID quality log: $!";
                my $fqual;
                while (! eof $qh) {
                    my $line = <$qh>;
                    if ($line =~ /Fine_Consistency=\s+(\d+(?:\.\d+)?)%/) {
                        $fqual = $1;
                    }
                }
                if (! $fqual) {
                    print "$genomeID quality check failed.\n";
                    $stats->Add(archaeaFail => 1);
                } else {
                    print "$genomeID is $fqual consistent.\n";
                    if ($fqual && $fqual >= 85) {
                        $good{$genomeID} = 1;
                        $stats->Add(archaeaGood => 1);
                    } else {
                        $stats->Add(archaeaBad => 1);
                    }
                }
            }
        }
    }
    # Now we run through the couplets and output any couplet for which $good{$genomeID} == 1.
    for my $couplet (@$couplets) {
        my ($genomeID, $line) = @$couplet;
        if ($good{$genomeID} && $good{$genomeID} == 1) {
            P3Utils::print_cols($line, oh => $oh);
        }
    }
}