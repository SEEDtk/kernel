use strict;
use FIG_Config;
use File::Copy::Recursive;
use GEO;
use EvalCon;
use Stats;
use P3DataAPI;

my ($inFile, $inDir) = @ARGV;
my $stats = Stats->new();
print STDERR "Input file is $inFile.  Input directory is $inDir.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
my $p3 = P3DataAPI->new();
my %geoOptions = (rolesHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, detail => 0, logH => \*STDERR);
# This is the original file.  We will create a new one with two additional columns.
open(my $ih, '<', $inFile) || die "Could not open input: $!";
# Skip the header and print the new header.
my $line = <$ih>;
print "Sample\tgood_meta\tbad_meta\tgood_patric\tbad_patric\tgood_bam\tbad_bam\n";
while (! eof $ih) {
    $line = <$ih>;
    chomp $line;
    my ($sample, $goodMeta, $badMeta, $goodPat, $badPat) = split /\t/, $line;
    $stats->Add(lineIn => 1);
    $stats->Add(goodMeta => $goodMeta);
    $stats->Add(allMeta => ($goodMeta + $badMeta));
    $stats->Add(goodPat => $goodPat);
    $stats->Add(allPat => ($goodPat + $badPat));
    # Get all the GTOs in the sample directory.
    my ($goodBam, $badBam) = (0, 0);
    my $sampDir = "$inDir/metabat-$sample";
    if (! -d $sampDir) {
        print STDERR "WARNING: $sampDir not found.\n";
        $stats->Add(missingDir => 1);
    } else {
        print STDERR "Processing $sampDir.\n";
        opendir(my $dh, $sampDir) || die "Could not open $sampDir: $!";
        my @gtoFiles = map { "$sampDir/$_" } grep { $_ =~ /\.gto$/ } readdir $dh;
        print STDERR scalar(@gtoFiles) . " genomes found.\n";
        # Create the GEOs.
        my $geoHash = GEO->CreateFromGtoFiles(\@gtoFiles, %geoOptions);
        my $geoCount = scalar keys %$geoHash;
        print STDERR "$geoCount genomes loaded.\n";
        $stats->Add(geosRead => $geoCount);
        for my $genome (sort keys %$geoHash) {
            my $geo = $geoHash->{$genome};
            if ($geo->is_good()) {
                $goodBam++;
            } else {
                $badBam++;
            }
        }
    }
    $stats->Add(goodBam => $goodBam);
    $stats->Add(allBam => ($goodBam + $badBam));
    print join("\t", $sample, $goodMeta, $badMeta, $goodPat, $badPat, $goodBam, $badBam) . "\n";
}
print STDERR "All done:\n" . $stats->Show();
