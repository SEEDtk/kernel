use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

$| = 1;
my $opt = ScriptUtils::Opts('dir');
my ($dir) = @ARGV;
opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
my @packages = grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
for my $package (@packages) {
    open(my $ih, '<', "$dir/$package/data.tbl") || die "Could not open data file for $package: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /\t$/) {
            print "Anomaly in $package.\n";
        }
    }
}
