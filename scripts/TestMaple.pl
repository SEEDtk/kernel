use strict;
use FIG_Config;
use File::Copy::Recursive;

if (-d "FakeCore") {
    File::Copy::Recursive::pathempty("FakeCore");
} else {
    File::Copy::Recursive::pathmk("FakeCore");
}
opendir(my $dh, "/vol/core-seed/FIGdisk/FIG/Data/Organisms") || die "Could not open FakeCore: $!";
my @orgs = grep { $_ =~ /^\d+\.\d+/ && -s "/vol/core-seed/FIGdisk/FIG/Data/Organisms/$_/assigned_functions" } readdir $dh;
print scalar(@orgs) . " genomes found.\n";
closedir $dh;
for my $org (@orgs) {
    my $orgDir = "/vol/core-seed/FIGdisk/FIG/Data/Organisms/$org";
    for my $file (map { "Features/$_/deleted.features" } qw(peg rna)) {
        if (-s "$orgDir/$file") {
            print "Copying $file.\n";
            File::Copy::Recursive::fcopy("$orgDir/$file", "FakeCore/$org/$file") || die "Could not copy $org $file: $!";
        }
        print "Copying $org.\n";
        File::Copy::Recursive::fcopy("$orgDir/assigned_functions", "FakeCore/$org/assigned_functions") || die "Could not copy $org functions: $!";
        File::Copy::Recursive::fcopy("$orgDir/GENOME", "FakeCore/$org/GENOME") || die "Could not copy $org genome name: $!";
    }
}
print "All done.\n";
