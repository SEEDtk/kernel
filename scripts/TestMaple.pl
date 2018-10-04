use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

$| = 1;
my %binCounts;
my $stats = Stats->new();
my $binDir = "$FIG_Config::data/Bins_HMP";
opendir(my $dh, $binDir) || die "Could not open binning directory: $!";
my @samples = grep { -s "$binDir/$_/contigs.fasta" } readdir $dh;
my $nSamples = scalar @samples;
print "$nSamples directories found.\n";
closedir $dh;
for my $sample (@samples) {
    $stats->Add(sampleIn => 1);
    my $sampleDir = "$binDir/$sample";
    if (-s "$sampleDir/site.tbl") {
        $stats->Add(siteFound => 1);
    } else {
        print "Setting site for $sample.\n";
        open(my $oh, ">$sampleDir/site.tbl") || die "Could not open site.tbl: $!";
        print $oh join("\t", 'RAE', 'unspecified', 'Unspecified') . "\n";
        $stats->Add(siteAdded => 1);
    }
}
print "All done.\n" . $stats->Show();