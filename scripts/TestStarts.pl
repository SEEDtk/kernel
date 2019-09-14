use strict;
use FIG_Config;
use IPC::Run3;
no warnings;

## Positional parameter is directory of GTO files.
my ($gtoDir) = @ARGV;
opendir(my $dh, $gtoDir) || die "Could not open input directory: $!";
my @gFiles = map { "$gtoDir/$_" } grep { $_ =~ /^\d+\.\d+\.gto$/ } readdir $dh;
closedir $dh;
for my $gto (@gFiles) {
    print STDERR "Processing $gto.\n";
    my $gName = "test";
    if ($gto =~ /(\d+\.\d+)/) {
        $gName = $1;
    }
    my $phase1in = "$FIG_Config::data/$gName.test1.tbl";
    my $phase1out = "$FIG_Config::data/$gName.pred1.tbl";
    my $phase2in =  "$FIG_Config::data/$gName.test2.tbl";
    my $phase2out = "$FIG_Config::data/$gName.pred2.tbl";
    run3 "genome.contigs test -u 45 -d 44 -f -v $gto", \undef, $phase1in;
    run3 "dl4j.run predict --meta location,codon,expect $FIG_Config::data/Coding/ComboTest", $phase1in, $phase1out;
    run3 "genome.starts test -v -r $FIG_Config::data/Global/roles.ser $gto $phase1out", \undef, $phase2in;
    run3 "dl4j.run predict --meta location,expect,roles $FIG_Config::data/Coding/FancyStartTest", $phase2in, $phase2out;
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
    print "$total predictions for $gName.\n";
    print "$hit hits. $miss misses.  $badMiss bad misses. sensitivity = " . ($hit / ($hit + $miss)) . "\n";
    print "$falseHit false hits.  specificity = " . ($trueN / ($trueN + $falseHit)) . "\n";
    print "Accuracy = " . (($trueN + $hit) / $total) . "\n";
}
