use strict;
use FIG_Config;

$| = 1;
my ($dir) = @ARGV;
opendir(my $dh, $dir) || die "Could not open $dir: $!";
my @bins = grep { -s "$dir/$_/run.log" } readdir $dh;
closedir $dh;
my ($count, $tot) = (0, scalar @bins);
print STDERR "$tot bins found.\n";
print "Bin\tsample\tfootprint\n";
for my $bin (@bins) {
    $count++;
    print STDERR "Processing $bin ($count of $tot).\n";
    open(my $ih, '<', "$dir/$bin/run.log") || die "Could not open log for $bin: $!";
    my ($size, $done, $footprint) = ('', 0, 0);
    while (! $done) {
        my $line = <$ih>;
        if (! defined $line || $line =~ /^SPAdes log can be found here/) {
            $done = 1;
        } elsif ($line =~ /\s+\d+G \/ (\d+)G\s+INFO/) {
            $footprint = $1 if $footprint < $1;
        } elsif ($line =~ /^Sample size is (\d+) gigabytes/) {
            $size = $1;
        }
    }
    if ($size) {
        print "$bin\t$size\t$footprint\n";
    }
}
