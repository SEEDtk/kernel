use strict;
use FIG_Config;
use File::Copy::Recursive;
use Bin;
use SeedUtils;

opendir(my $dh, 'Bins_HMP') || die "Could not open binning directory: $!";
my @samples = grep { -s "Bins_HMP/$_/contigs.fasta" } readdir $dh;
closedir $dh;
for my $sample (@samples) {
    print STDERR "Processing $sample.\n";
    my $subDir = "Bins_HMP/$sample";
    my $target = "Bins_Test/$sample";
    my @files = qw(contigs.fasta site.tbl output.contigs2read.txt);
    for my $file (@files) {
        File::Copy::Recursive::fcopy("$subDir/$file", "$target/$file");
    }
}