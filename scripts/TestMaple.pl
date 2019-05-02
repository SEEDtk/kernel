use strict;
use FIG_Config;
use GeoGroup;
use P3Utils;
use Stats;

$| = 1;
my $rootDir = "$FIG_Config::data/Italians";
opendir(my $dh, $rootDir) || die "Could not open italian directory: $!";
my @dirs = grep { -d "$rootDir/$_" } readdir $dh;
closedir $dh;
my $stats = Stats->new();
P3Utils::print_cols(['PATRIC ID', 'Good Seed', 'Good']);
for my $dir (@dirs) {
    print STDERR "Processing directory $dir.\n";
    $stats->Add(dirIn => 1);
    my $geoGroup = GeoGroup->new({stats => $stats, logH => \*STDERR, detail => 0}, "$rootDir/$dir");
    my $geoList = $geoGroup->geoList;
    for my $geo (@$geoList) {
        my ($goodSeed, $good) = ('', '');
        if ($geo->good_seed) {
            $goodSeed = 'Y';
            $stats->Add(goodSeed => 1);
        } else {
            $stats->Add(badSeed => 1);
        }
        if ($geo->is_good) {
            $good = 'Y';
            $stats->Add(goodGenome => 1);
        } else {
            $stats->Add(badGenome => 1);
        }
        P3Utils::print_cols([$geo->id, $goodSeed, $good]);
        $stats->Add(genomeIn => 1);
    }
}
print STDERR "All done.\n" . $stats->Show();
