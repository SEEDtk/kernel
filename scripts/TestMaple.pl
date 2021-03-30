use strict;
use FIG_Config;
use File::Copy::Recursive;

opendir(my $dh, "Bins_HMP") || die "Could not open Bins_HMP: $!";
my @dirs = grep { -s "Bins_HMP/$_/Eval/index.tbl" } readdir $dh;
print "sample_id\tgood\tbad\n";
for my $dir (@dirs) {
    print STDERR "Processing $dir.\n";
    my $subDir = "Bins_HMP/$dir";
    if (open(my $ih, '<', "$subDir/site.tbl")) {
        my $line = <$ih>;
        close $ih; undef $ih;
        if ($line =~ /HMP\s+stool/i) {
            open($ih, '<', "$subDir/Eval/index.tbl") || die "Could not open index.tbl for $dir: $!";
            $line = <$ih>;
            my ($good, $bad) = (0, 0);
            while (! eof $ih) {
                $line = <$ih>;
                if ($line =~ /1$/) {
                    $good++;
                } else {
                    $bad++;
                }
            }
            print "$dir\t$good\t$bad\n";
        }
    }
}
