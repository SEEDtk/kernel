use strict;
use FIG_Config;
use ScriptUtils;
use SeedUtils;
use Data::Dumper;
use Proc::ParallelLoop;

my @inputs = qw(10 12 15 18 20 22 25 28 30 32 35 38 40 42 45 48 50 52 55 58
                60 62 65 68 70 72 75 78 80 82 85 88 90 92 95 98 100);
my @scores;
print "Processing...\n";
Proc::ParallelLoop::pareach(\@inputs, \&process);
print scalar(@scores) . " scores found.\n";
for my $score (@scores) {
    print "$score ";
}
print "\n";


sub process {
    my ($a) = @_;
    if ($a % 3 == 0) {
        push @scores, $a
    }
}
