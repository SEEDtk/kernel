use strict;
use FIG_Config;
use IPC::Run3;
no warnings;

## Positional parameter is name of GTO file.
my ($gto) = @ARGV;
my $gName = "test";
if ($gto =~ /(\d+\.\d+)/) {
	$gName = $1;
}
my $phase1in = "$FIG_Config::data/$gName.test1.tbl";
my $phase1out = "$FIG_Config::data/$gName.pred1.tbl";
my $phase2in =  "$FIG_Config::data/$gName.test2.tbl";
my $phase2out = "$FIG_Config::data/$gName.pred2.tbl";
run3 ["genome.contigs", qw(test -u 45 -d 44 -f -v), $gto], \undef, $phase1in;
run3 ["dl4j.run", qw(predict --meta location,codon,expect), "$FIG_Config::data/Coding/ComboTest"], $phase1in, $phase1out;
run3 ["genome.starts", qw(test -v -r), "$FIG_Config::data/Global/roles.ser", $gto, $phase1out], \undef, $phase2in;
run3 ["dl4j.run", qw(predict --meta location,expect,roles), "$FIG_Config::data/Coding/FancyStartTest"], $phase2in, $phase2out;
print STDERR "Analyzing predictions in $phase2out.\n";
open(my $ih, '<', $phase2out) || die "Could not open predictions file: $!";
my $line = <$ih>;
my ($hit, $miss, $trueN, $falseHit, $badMiss, $total) = (0, 0, 0, 0, 0, 0);
while (! eof $ih) {
	my $line = <$ih>;
	my ($loc, $expect, $roles, $predicted, $conf) = split /\t/, $line;
	if ($expect eq 'start') {
		if ($predicted eq 'start') {
			$hit++;
		} else {
			$miss++;
			if ($roles) {
				$badMiss++;
			}
		}
	} else {
		if ($predicted eq 'other') {
			$trueN++;
		} else {
			$falseHit++;
		}
	}
	$total++;
}
close $ih;
print STDERR "$total predictions.\n";
print STDERR "$hit hits. $miss misses.  $badMiss bad misses. sensitivity = " . ($hit / ($hit + $miss)) . "\n";
print STDERR "$falseHit false hits.  specificity = " . ($trueN / ($trueN + $falseHit)) . "\n";
print STDERR "Accuracy = " . (($trueN + $hit) / $total) . "\n";
