use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use Shrub;
use Contigs;
use SeedUtils;

opendir(my $dh, "Bins_HMP") || die "Could not open bin directory: $!";
my @bins = grep { -s "Bins_HMP/$_/bins.found.tbl" } readdir $dh;
closedir $dh;
my $found = scalar @bins;
print "$found bins found.\n";
my $count = 0;
open(my $oh, '>', "samples.fna") || die "Could not open samples output: $!";
open(my $ph, '>', "samples.faa") || die "Could not open samples output: $!";
for my $bin (@bins) {
    $count++;
    print "Reading contigs for $bin ($count of $found).\n";
    my $contigs = Contigs->new("Bins_HMP/$bin/contigs.fasta", genomeID => $bin);
    open(my $ih, '<', "Bins_HMP/$bin/bins.found.tbl") || die "Could not open $bin input: $!";
    my $subseq = 1;
    while (my $line = <$ih>) {
        chomp $line;
        my ($contig, $start, $dir, $len) = split /\t/, $line;
        if ($dir eq '-') {
            $start -= ($len - 1);
        }
        print "$subseq is $contig, $start $dir $len.\n";
        my $dna = $contigs->dna([$contig, $start, $dir, $len]);
        my $seq = $contigs->xlate([$contig, $start, $dir, $len]);
        print $ph ">$bin.$contig $subseq\n$seq\n";
        print $oh ">$bin.$contig $subseq\n$dna\n";
        $subseq++;
    }
    print "$subseq sequences output for $bin.\n";
}
close $oh;