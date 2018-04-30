use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use GPUtils;

$| = 1;
my $stats = Stats->new();
my $gHash = GPUtils::get_all('ModPackages');
for my $genome (sort keys %$gHash) {
    print STDERR "Processing $genome.\n";
    my $gto = GPUtils::gto_of($gHash, $genome);
    $stats->Add(gtoRead => 1);
    my $flag = (GPUtils::good_seed($gto) ? 1 : 0);
    if ($flag) {
        $stats->Add(gtoGood => 1);
    }
    print "$genome\t$flag\n";
}
print STDERR "All done.\n" . $stats->Show();