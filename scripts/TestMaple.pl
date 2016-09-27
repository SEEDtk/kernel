use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

$| = 1;
my $opt = ScriptUtils::Opts('gfile');
my ($qfile) = @ARGV;
open(my $ih, "<$qfile") || die "Could not open $qfile: $!";
my $mode = "";
my $stats = Stats->new();
while (! eof $ih) {
    my $line = <$ih>;
    my $count = $stats->Add(lines => 1);
    print "$count lines read.\n" if ($count % 10000 == 0);
    if ($mode eq '@') {
        $stats->Add(reads => 1);
        chomp $line;
        my $n = length($line);
        my $ns = ($line =~ tr/ACGT//);
        my $xs = ($line =~ tr/N//);
        if ($ns == $n) {
            $stats->Add(goodLine => 1);
        } elsif ($ns == 0) {
            $stats->Add(badLine => 1);
        }
        if ($xs == $n) {
            $stats->Add(nLine => 1);
        }
        $stats->Add(nChar => $xs);
        $stats->Add(goodChar => $ns);
        $stats->Add(badChar => ($n - $ns));
        $mode = "";
    } else {
        if (substr($line,0,1) eq '@') {
            $mode = '@';
        } else {
            $mode = "";
        }
    }
}
print "Counts:\n" . $stats->Show();
