use strict;
use FIG_Config;
use File::Copy::Recursive;
use Bin;
use SeedUtils;

opendir(my $dh, 'Bins_HMP') || die "Could not open binning directory: $!";
my @samples = grep { -s "Bins_HMP/$_/bins.rast.json" } readdir $dh;
closedir $dh;
for my $sample (@samples) {
    print STDERR "Processing $sample.\n";
    my $error = 0;
    my $subDir = "Bins_HMP/$sample";
    my $binList = Bin::ReadBins("$subDir/bins.rast.json");
    for (my $idx = 1; $idx <= scalar @$binList && ! $error; $idx++) {
        my $bin = $binList->[$idx - 1];
        my $gtoFile = "$subDir/bin$idx.gto";
        if (! -s $gtoFile) {
            $error++;
            print STDERR "$gtoFile is missing for bin $idx.\n";
        } else {
            my $gto = SeedUtils::read_encoded_object($gtoFile);
            my $gtoName = $gto->{scientific_name};
            $gtoName =~ s/\s+cleaned//;
            my $binName = $bin->name || "(unknown)";
            if ($gtoName ne $binName) {
                print STDERR "Bin $idx should be $binName but is $gtoName.\n";
                $error++;
            }
        }
    }
    if ($error) {
        print "$sample\n";
    }
}