use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

$| = 1;
my %binCounts;
my $stats = Stats->new();
open(my $bh, "<$FIG_Config::data/binList.tbl") || die "Could not open binList.tbl: $!";
# File is sampleID, good, bad.
while (! eof $bh) {
    my $line = <$bh>;
    if ($line =~ /(\S+)\t(\d+)\t(\d+)/) {
        $binCounts{$1} = [$2, $3];
        $stats->Add(countsIn => 1);
    }
}
close $bh;
my $binDir = "$FIG_Config::data/Bins_HMP";
opendir(my $dh, $binDir) || die "Could not open binning directory: $!";
my @samples = grep { -s "$binDir/$_/contigs.fasta" } readdir $dh;
my $nSamples = scalar @samples;
print STDERR "$nSamples directories found.\n";
print join("\t", qw(Sample GoodBins BadBins Coverage ContigLength ReadLength)) . "\n";
my $kSamples = 0;
closedir $dh;
for my $sample (@samples) {
    $kSamples++;
    $stats->Add(sampleIn => 1);
    print STDERR "Processing $sample ($kSamples of $nSamples).\n";
    open(my $fh, "<$binDir/$sample/contigs.fasta") || die "Could not open contigs.fasta for $sample: $!";
    my $len = 0;
    my $cov = 0;
    while (! eof $fh) {
        my $line = <$fh>;
        $stats->Add(fastaLineIn => 1);
        if ($line =~ /^>NODE_\d+_length_(\d+)_cov_(\d+)/) {
            $len += $1;
            $cov += ($1 * $2);
            $stats->Add(fastaLineUsed => 1);
        }
    }
    if ($len == 0) {
        print STDERR "Sample $sample was not well-formed.\n";
        $stats->Add(sampleRejected => 1);
    } else {
        my $counts = $binCounts{$sample} // [0, 0];
        my ($good, $bad) = @$counts;
        my $realCov = $cov / $len;
        print join("\t", $sample, $good, $bad, $realCov, $len, int($cov)) . "\n";
    }
}
print STDERR "All done.\n" . $stats->Show();