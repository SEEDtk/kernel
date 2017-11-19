use strict;
use FIG_Config;
use ScriptUtils;
use Shrub;
use Stats;
use SeedUtils;

my ($dir) = @ARGV;
opendir(my $dh, $dir);
my @subs = grep { -f "$dir/$_/contigs.fasta" } readdir $dh;
for my $sub (@subs) {
    mkdir "$dir/$sub/Assembly";
    print "$sub fixed.\n";
}