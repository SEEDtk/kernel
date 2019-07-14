use strict;
use FIG_Config;
use File::Copy::Recursive;
use Bin;
use SeedUtils;

opendir(my $dh, 'Bins_HMP') || die "Could not open binning directory: $!";
my @samples = grep { -s "Bins_HMP/$_/Eval/index.tbl" } readdir $dh;
closedir $dh;
for my $sample (@samples) {
    print STDERR "Processing $sample.\n";
    my $subDir = "Bins_HMP/$sample";
    open(my $ih, '<', "$subDir/Eval/index.tbl") || die "Could not open index.tbl: $!";
    my $line = <$ih>;
    my %good;
    while (! eof $ih) {
        $line = <$ih>;
        if ($line =~ /(\d+\.\d+).+1$/) {
            $good{$1} = 1;
        }
    }
    my $idx = 1;
    while (-s "$subDir/bin$idx.gto") {
        my $gto = SeedUtils::read_encoded_object("$subDir/bin$idx.gto");
        my $id = $gto->{id};
        if ($good{$id}) {
            print "Copying bin $idx: $id\n";
            File::Copy::Recursive("$subDir/bin$idx.gto", "GoodBins/$id.gto");
        }
    }
}