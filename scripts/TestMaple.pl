use strict;
use FIG_Config;
use File::Copy::Recursive;

my ($inDir) = @ARGV;
opendir(my $dh, $inDir) || die "Could not open input directory: $!";
my @files = grep { $_ =~ /\.fna$/ } readdir $dh;
closedir $dh;
my $num = 1;
my $fileCount = 0;
my $lineCount = 0;
my $changeCount = 0;
for my $file (@files) {
    $fileCount++;
    open(my $ih, '<', "$inDir/$file") || die "Could not open $file: $!";
    open(my $oh, '>', "$inDir/temp.fna") || die "Could not open output file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        $lineCount++;
        if ($line =~ /^>REP_(\d+\.\d+_)\d+(_.+)/) {
            my $label = ">REP_$1$num$2";
            print $oh "$label\n";
            $changeCount++;
            $num++;
        } else {
            print $oh $line;
        }
    }
    close $oh; close $ih;
    File::Copy::Recursive::fmove("$inDir/temp.fna", "$inDir/$file") || die "Could not rename temp file: $!";
}
print "$fileCount files, $lineCount lines, $changeCount changed.\n";
