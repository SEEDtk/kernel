use strict;
use FIG_Config;
use SRAlib;
use Stats;

$| = 1;
open(my $ih, '<', 'sites.tbl') || die "Could not open input file: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($sample, $project, $site) = split /\t/, $line;
    open(my $oh, '>', "Bins_HMP/$sample/site.tbl") || die "Could not open site.tbl for $sample: $!";
    my $siteTitle = join(' ', map { ucfirst $_ } split /_/, $site);
    print join("\t", $project, $site, $siteTitle) . "\n";
    close $oh;
}