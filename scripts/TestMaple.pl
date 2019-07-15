use strict;
use FIG_Config;
use File::Copy::Recursive;
use Bin;
use SeedUtils;

my @files = qw(contigs.fasta output.contigs2reads.txt site.tbl);
opendir(my $dh, 'Bins_Test') || die "Could not open binning directory: $!";
my @samples = grep { -s "Bins_Test/$_/contigs.fasta" } readdir $dh;
closedir $dh;
my $count = 0;
for my $sample (@samples) {
    print STDERR "Processing $sample.\n";
    my $subDir = "Bins_Test/$sample";
    my $oldDir = "Bins_HMP/$sample";
    File::Copy::Recursive::pathempty($subDir) || die "Could not erase $subDir: $!";
    for my $file (@files) {
        File::Copy::Recursive::fcopy("$oldDir/$file", "$subDir/$file") || die "Could not copy $file into $subDir: $!";
    }
}