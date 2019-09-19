use strict;
use FIG_Config;
use SeedUtils;
use File::Copy::Recursive;

my @files = qw(GENOME TAXONOMY GENETIC_CODE TAXONOMY_ID);
my $orgDir = '/vol/core-seed/FIGdisk/FIG/Data/Organisms';
opendir(my $dh, $orgDir) || die "Could not open coreSEED directory: $!";
my @genomes = grep { -d "$orgDir/$_" && $_ =~ /^(\d+\.\d+)$/ } readdir $dh;
closedir $dh;
for my $genome (@genomes) {
    print "Processing $genome.\n";
    File::Copy::Recursive::pathmk("FakeCoreOrgs/$genome") || die "Could not create $genome directory: $!";
    for my $file (@files) {
        if (-f "$orgDir/$genome/$file") {
            File::Copy::Recursive::fcopy("$orgDir/$genome/$file", "FakeCoreOrgs/$genome") || die "Could not copy $genome $file: $!";
        }
    }
}
