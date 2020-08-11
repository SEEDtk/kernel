use strict;
use FIG_Config;
use File::Copy::Recursive;

opendir(my $dh, "FakeCore") || die "Could not open FakeCore: $!";
my @orgs = grep { $_ =~ /^\d+\.\d+/ } readdir $dh;
closedir $dh;
for my $org (@orgs) {
    my $orgDir = "/vol/core-seed/FIGdisk/FIG/Data/Organisms/$org";
    for my $file (map { "Features/$_/deleted.features" } qw(peg rna)) {
        if (-s "$orgDir/$file") {
            print "Copying $file.\n";
            File::Copy::Recursive::fcopy("$orgDir/$file", "FakeCore/$file") || die "Could not copy $org $file: $!";
        }
    }
}
print "All done.\n";
