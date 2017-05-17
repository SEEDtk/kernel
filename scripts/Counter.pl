use Shrub;

my $shrub = Shrub->new();
my $count = $shrub->GetCount('Protein', '', []);
print "$count proteins.\n";
$count = $shrub->GetCount('Genome', '', []);
print "$count genomes.\n";
$count = $shrub->GetCount('Subsystem', '', []);
print "$count subsystems.\n";
