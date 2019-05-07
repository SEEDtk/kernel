use strict;
use FIG_Config;
use GEO;
use P3Utils;
use Stats;

$| = 1;
my %loaded;
open(my $ih, '<', 'temp2.tbl') || die "Could not open temp2.tbl: $!";
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /(\d+\.\d+)/) {
        $loaded{$1} = 1;
    }
}
print scalar(keys %loaded) . " already loaded.\n";
close $ih;
open(my $oh, '>>', 'temp2.tbl') || die "Could not open temp2.tbl: $!";
my $rootDir = "$FIG_Config::data/Italians";
opendir(my $dh, $rootDir) || die "Could not open italian directory: $!";
my @dirs = grep { substr($_,0,1) ne '.' && -d "$rootDir/$_" } readdir $dh;
closedir $dh; undef $dh;
my $stats = Stats->new();
for my $dir (@dirs) {
    print "Processing directory $dir.\n";
    $stats->Add(dirIn => 1);
    my $subDir = "$rootDir/$dir";
    opendir($dh, $subDir) || die "Could not open $subDir: $!";
    my @files = grep { $_ =~ /^(\d+\.\d+).gto$/ && ! $loaded{$1} } readdir $dh;
    if (! @files) {
        print "Directory already loaded.\n";
    } else {
        my $gHash = GEO->CreateFromGtoFiles([ map { "$subDir/$_" } @files ], stats => $stats, logH => \*STDOUT, detail => 0);
        for my $genome (keys %$gHash) {
            my $geo = $gHash->{$genome};
            my ($goodSeed, $good) = ('', '');
            if ($geo->good_seed) {
                $goodSeed = 'Y';
                $stats->Add(goodSeed => 1);
            } else {
                $stats->Add(badSeed => 1);
            }
            if ($geo->is_good) {
                $good = 'Y';
                $stats->Add(goodGenome => 1);
            } else {
                $stats->Add(badGenome => 1);
            }
            P3Utils::print_cols([$genome, $goodSeed, $good], oh => $oh);
            $stats->Add(genomeIn => 1);
        }
    }
}
print "All done.\n" . $stats->Show();
