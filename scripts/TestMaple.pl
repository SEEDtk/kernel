use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

my $opt = ScriptUtils::Opts('gfile');
my ($qfile) = @ARGV;
open(my $ih, "<$qfile") || die "Could not open $qfile: $!";
my $mode = "";
my $stats = Stats->new();
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(lines => 1);
    if ($mode eq '@') {
        $stats->Add(reads => 1);
        chomp $line;
        my $n = length($line);
        my $ns = ($line =~ tr/ACGT//);
        if ($ns == $n) {
            $stats->Add(goodLine => 1);
        } elsif ($ns == 0) {
            $stats->Add(badLine => 1);
        }
        $stats->Add(goodChar => $ns);
        $stats->Add(badChar => ($n - $ns));
    } else {
        if (substr($line,0,1) eq '@') {
            $mode = '@';
        } else {
            $mode = "";
        }
    }
}
print "Counts:\n" . $stats->Show();
