use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use FastA;
use Number::Format;

my $f = new Number::Format(decimal_digits => 0, thousands_sep => ',');
my %binCounts; # #bins -> [count,min,max,tot]
my $binDir = "$FIG_Config::data/Bins_HMP";
opendir(my $dh, $binDir) || die "Could not open bins directory: $!";
my @samples = sort grep { -s "$binDir/$_/bins.report.txt" } readdir $dh;
for my $sample (@samples) {
    open(my $ih, "<$binDir/$sample/bins.report.txt") || die "Could not open bin report: $!";
    my $binCount = 0;
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^BIN \d/) {
            $binCount++;
        }
    }
    close $ih; undef $ih;
    my $bp = 0;
    if (open($ih, "<$binDir/$sample/contigs.fasta")) {
        my $reader = FastA->new($ih);
        while ($reader->next) {
            my $seq = $reader->left;
            $bp += length $seq;
        }
    }
    print STDERR "$sample has $binCount bins and $bp base pairs.\n";
    my $countData = $binCounts{$binCount};
    if (! $countData) {
        $binCounts{$binCount} = [1, $bp, $bp, $bp];
    } else {
        $countData->[0]++;
        $countData->[1] = $bp if $bp < $countData->[1];
        $countData->[2] = $bp if $bp > $countData->[2];
        $countData->[3] += $bp;
    }
}
print STDERR "Writing results.\n";
print "#bins\tnum\tmin\tmax\tmean\n";
for my $count (sort { $a <=> $b } keys %binCounts) {
    my ($num, $min, $max, $tot) = @{$binCounts{$count}};
    my @row = map { $f->format_number($_) } ($count, $num, $min, $max, $tot/$num);
    print join("\t", @row) . "\n";
}
