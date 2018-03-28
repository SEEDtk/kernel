use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use Shrub;

my $shrub = Shrub->new();
my @taxa = qw(1100 1120941 1236494 411483);
for my $taxon (@taxa) {
    my @taxList = $shrub->taxonomy_of($taxon);
    print "$taxon is: " . join(", ", @taxList) . "\n";
}