use strict;
use FIG_Config;
use FastA;
use SeedUtils;
use Bin;

$| = 1;
my $rootDir = "$FIG_Config::data/Bins_HMP";
opendir(my $dh, $rootDir) || die "Could not open binning directory: $!";
my @dirs = grep { -s "$rootDir/$_/bins.rast.json" } readdir $dh;
closedir $dh;
print join("\t", qw(sample id name complete contam group multi)) . "\n";
for my $dir (@dirs) {
    # This will hold the taxon IDs of the bins with multiple reference genomes.
    my %multiTaxes;
    if (! open(my $ih, '<', "$rootDir/$dir/bins.rast.json")) {
        print STDERR "Open failed for bin file of $dir: $!\n";
    } else {
        print STDERR "Processing bin file for $dir.\n";
        my $binList = Bin::ReadBins($ih);
        for my $bin (@$binList) {
            my @refs = $bin->refGenomes;
            if (scalar(@refs) > 1) {
                $multiTaxes{$bin->taxonID} = 1;
            }
        }
    }
    # Now loop through the evaluation index.  We want bins that are not good, have high completeness, but also high
    # contamination.
    if (! open(my $kh, '<', "$rootDir/$dir/Eval/index.tbl")) {
        print STDERR "Open failed for quality file of $dir: $!\n";
    } else {
        my $needsGroups = 0;
        print STDERR "Processing quality file for $dir.\n";
        my $line = <$kh>;
        while (! eof $kh) {
            $line = <$kh>;
            $line =~ s/[\r\n]+$//;
            my @fields = split /\t/, $line;
            my ($id, $name, $complt, $contam, $group) = ($fields[1], $fields[2], $fields[10], $fields[11], $fields[12]);
            if ($complt >= 90 && $contam > 10) {
                my ($tax) = split /\./, $id;
                my $multiFlag = ($multiTaxes{$tax} ? 'multi' : '');
                print join("\t", $dir, $id, $name, $complt, $contam, $group, $multiFlag) . "\n";
            }
        }
    }
}