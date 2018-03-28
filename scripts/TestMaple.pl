use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use Shrub;

my $shrub = Shrub->new();
my @taxa = qw(1100 1120941 1236494 411483);
for my $taxon (@taxa) {
    my $done;
    while (! $done) {
        my ($taxData) = $shrub->GetAll('IsInTaxonomicGroup TaxonomicGrouping', 'IsInTaxonomicGroup(from-link) = ?', [$taxon],
            'TaxonomicGrouping(id) TaxonomicGrouping(name) TaxonomicGrouping(domain)');
        if (! $taxData) {
            print "No parent found for $taxon.\n";
            $done = 1;
        } else {
            my ($id, $name, $domainFlag) = @$taxData;
            print "Parent of $taxon is $id: $name.\n";
            $done = $domainFlag;
        }
    }
    print "\n";
}