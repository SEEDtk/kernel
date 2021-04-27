use strict;
use FIG_Config;
use File::Copy::Recursive;

my ($dest) = @ARGV;
opendir(my $dh, "Bins_HMP") || die "Could not open Bins_HMP: $!";
my @dirs = grep { -s "Bins_HMP/$_/contigs.fasta" } readdir $dh;
for my $dir (@dirs) {
    if (! open(my $ih, '<', "Bins_HMP/$dir/site.tbl")) {
        print "$dir has no site file.\n";
    } else {
        my $line = <$ih>;
        close $ih;
        if ($line =~ /HMP\s+stool/) {
            File::Copy::Recursive::fcopy("Bins_HMP/$dir/contigs.fasta", "$dest/$dir.fasta")
                || die "Error copying $dir: $!";
            print "$dir copied.\n";
        } else {
            print "$dir is of type $line";
        }
    }
}
