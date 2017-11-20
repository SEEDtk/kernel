use strict;
use FIG_Config;
use ScriptUtils;
use Shrub;
use Stats;
use SeedUtils;

my ($dir) = @ARGV;
opendir(my $dh, $dir);
my @subs = grep { -f "$dir/$_/contigs.fasta" } readdir $dh;
for my $sub (@subs) {
    my $workDir = "$dir/$sub";
    opendir(my $dh, $workDir) || die "Could not open work directory: $!";
    my @files = grep { -f "$workDir/$_" } readdir $dh;
    print "Deleting intermediate files in $workDir.\n";
    my ($count, $total) = (0,0);
    for my $file (@files) {
        my $fullName = "$workDir/$file";
        $total++;
        unless ($fullName =~ /abundance/ || $fullName =~ /\.fastq$/ || $fullName =~ /\.fq/ || $file eq 'site.tbl' ||
            $file eq 'contigs.fasta' || $file eq 'output.contigs2reads.txt') {
            unlink $fullName;
            $count++;
        }
    }
    print "$count of $total files deleted.\n";
    print "$sub fixed.\n";
}