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
        for (my $i = 0; $i < $n; $i++) {
            my $c = substr($line,$i,1);
            $stats->Add("char-$c" => 1);
        }
    } else {
        if (substr($line,0,1) eq '@') {
            $mode = '@';
        } else {
            $mode = "";
        }
    }
}
print "Counts:\n" . $stats->Show();
