use strict;
use FIG_Config;
use SRAlib;
use Stats;

$| = 1;
my $stats = Stats->new();
my $sraLib = SRAlib->new(logH => \*STDERR, stats => $stats);
print join("\t", qw(sample project site title)) . "\n";
print STDERR "Reading sample directory.\n";
opendir(my $dh, 'Bins_HMP') || die "Could not open Bins_HMP: $!";
my @samples = grep { substr($_,0,1) ne '.' &&  -d "Bins_HMP/$_" } readdir $dh;
for my $sample (@samples) {
    my $siteFile = "Bins_HMP/$sample/site.tbl";
    if (-s $siteFile) {
        $stats->Add(siteKnown => 1);
        open(my $ih, '<', $siteFile) || die "Could not open $siteFile: $!";
        my $line = <$ih>;
        print "$sample\t$line";
    } else {
        $stats->Add(siteNeeded => 1);
        my ($project, $site) = $sraLib->compute_site($sample);
        my $siteTitle = join(' ', map { ucfirst $_ } split /_/, $site);
        my $line = "$project\t$site\t$siteTitle\n";
        open(my $oh, '>', $siteFile) || die "Could not open $siteFile: $!";
        print $oh $line;
        print "$sample\t$line";
    }
}
print STDERR "All done.\n" . $stats->Show();
