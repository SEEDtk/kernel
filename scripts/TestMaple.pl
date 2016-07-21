use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

my $stats = Stats->new();
my $opt = ScriptUtils::Opts('mainDir', ['test', 'trace without deleting']);
my ($mainDir) = @ARGV;
opendir(my $wh, $mainDir) || die "Could not open main directory: $!";
my @workDirs = grep { substr($_,0,1) ne '.' && -d "$mainDir/$_" } readdir $wh;
close $wh;
for my $work (@workDirs) {
    $stats->Add(totalDirs => 1);
    my $workDir = "$mainDir/$work";
    opendir(my $dh, $workDir) || die "Could not open work directory: $!";
    my @files = grep { -f "$workDir/$_" } readdir $dh;
    print "Deleting intermediate files in $workDir.\n";
    my ($count, $total) = (0,0);
    for my $file (@files) {
        my $fullName = "$workDir/$file";
        $stats->Add(totalFiles => 1);
        unless ($file =~ /_abundance_table\.tsv$/ || $file =~ /\.fastq$/ ||
                $file eq 'contigs.fasta' || $file eq 'output.contigs2reads.txt' || $file eq 'run.log') {
            if ($opt->test) {
                print "Delete $fullName\n";
            } else {
                unlink $fullName;
            }
            $stats->Add(deletedFiles => 1);
        }
    }
}
print "All done.\n" . $stats->Show();
