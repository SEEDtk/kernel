use strict;
use FIG_Config;
use SRAlib;
use Stats;

$| = 1;
my $stats = Stats->new();
my $sraLib = SRAlib->new(logH => \*STDERR, stats => $stats);
print join("\t", qw(sample key value)) . "\n";
print STDERR "Reading sample directory.\n";
opendir(my $dh, 'Bins_HMP') || die "Could not open Bins_HMP: $!";
my @samples = grep { $_ =~ /^SRS\d+$/ && -d "Bins_HMP/$_" } readdir $dh;
for my $sample (@samples) {
    if (-s "Bins_HMP/$sample/site.tbl") {
        $stats->Add(siteKnown => 1);
    } else {
        $stats->Add(siteNeeded => 1);
        my $metaH = $sraLib->get_meta($sample);
        for my $key (sort keys %$metaH) {
            print join("\t", $sample, $key, $metaH->{$key}) . "\n";
        }
    }
}
print "All done.\n" . $stats->Show();
