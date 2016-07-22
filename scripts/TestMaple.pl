use strict;
use FIG_Config;
use ScriptUtils;
use File::Copy;

my $opt = ScriptUtils::Opts('sampleDir');
my ($sampleDir) = @ARGV;
my $inputDir = "$FIG_Config::data/InputsSamples";
my @subs = qw(HMP MH);
for my $sub (@subs) {
    my $subDir = "$inputDir/$sub";
    opendir(my $wh, $subDir) || die "Could not open $sub directory: $!";
    my @dirs = grep { substr($_,0,1) ne '.' && -d "$subDir/$_" } readdir $wh;
    for my $dir (@dirs) {
        copy("$subDir/$dir/site.tbl", "$sampleDir/$dir/site.tbl") || die "Copy failed for $dir: $!";
    }
}
