use strict;
use FIG_Config;
use GeoGroup;
use Stats;

$| = 1;
my $stats = Stats->new();
my %options = (stats => $stats, detail => 0, logH => \*STDERR);
P3Utils::print_cols(['score', 'Italian', 'Italian_name', 'Bin', 'Bin_name', 'Group']);
# Comparing Italian binning sets with ours.  Italians in IT_Bins; ours in Bins_IT.
opendir(my $th, 'IT_Bins') || die "Could not open italian directory: $!";
my @itDirs = grep { substr($_,0,1) ne '.' && -d "IT_Bins/$_" } readdir $th;
closedir $th;
for my $itDir (@itDirs) {
    print STDERR "Processing $itDir.\n";
    if (! -d "Bins_IT/$itDir") {
        print STDERR "No corresponding bin directory to $itDir.\n";
        $stats->Add(missingDir => 1);
    } else {
        $stats->Add(foundDir => 1);
        my $itGroup = GeoGroup->new(\%options, "IT_Bins/$itDir");
        my $binGroup = GeoGroup->new(\%options, "Bins_IT/$itDir");
        my ($pairs, $orphans1, $orphans2) = $itGroup->MapGroups($binGroup, 300);
        $stats->Add(itExtras => scalar @$orphans1);
        $stats->Add(binExtras => scalar @$orphans2);
        for my $pair (@$pairs) {
            my $itGeo = $itGroup->geo($pair->[0]);
            my $binGeo = $binGroup->geo($pair->[1]);
            my $output = 1;
            if (! $itGeo->is_good) {
                $stats->Add(badItGenome => 1);
                $output = 0;
            }
            if (! $binGeo->is_good) {
                $stats->Add(banBinGenome => 1);
                $output = 0;
            }
            if ($output) {
                P3Utils::print_cols([$pair->[2], $itGeo->id, $itGeo->name, $binGeo->id, $binGeo->name, $itDir]);
                $stats->Add(pairFound => 1);
            }
        }
    }
}
print STDERR "All done.\n" . $stats->Show();