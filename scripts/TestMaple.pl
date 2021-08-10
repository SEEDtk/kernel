use strict;
use FIG_Config;
use File::Copy::Recursive;

my $inDir = '/vol/patric3/QA/applications';
my $outDir = "$FIG_Config::data/p3Tests";
File::Copy::Recursive::pathempty($outDir) || die "Could not empty $outDir: $!";
opendir(my $dh, $inDir) || die "Could not open input directory: $!";
my @appdirs = grep { $_ =~ /^App-/ && -d "$inDir/$_/tests" } readdir $dh;
closedir $dh; undef $dh;
for my $appdir (@appdirs) {
    my $source = "$inDir/$appdir/tests";
    my $appName = substr $appdir, 4;
    print "Processing application $appName.\n";
    opendir(my $dh, $source) || die "Could not open $source directory: $!";
    my @jsons = grep { $_ =~ /\.json$/ } readdir $dh;
    close $dh; undef $dh;
    for my $json (@jsons) {
        File::Copy::Recursive::fcopy("$source/$json", "$outDir/$appName.$json") || die "Copy failed for $json of $appName: $!";
    }
}
print "All done.";
