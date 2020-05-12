use strict;
use FIG_Config;
use RoleParse;
use File::Copy::Recursive;

if (! -d "Profiles") {
    File::Copy::Recursive::pathmk("Profiles");
}
open(my $rh, '<', "Global/uniRoles.tbl") || die "Could not open uniRoles: $!";
my %checksum;
print "Reading roles.\n";
while (! eof $rh) {
    my ($rid, $check, $name) = split /\t/, <$rh>;
    $checksum{$check} = $rid;
}
close $rh;
print scalar(keys %checksum) . " universal roles found.\n";
print "Reading census.\n";
my $dir = "/disks/ssd/olson/profiles/bob1/";
open(my $ih, '<', "/homes/gdpusch/Projects/Profiles/cluster_census.subsystem_based.tab") || die "Could not open census: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($subDir, $name, $pMembers, $members, $percent, $function) = split /\t/, $line;
    if ($members >= 10 && $pMembers / $members >= 0.9) {
        my $check = RoleParse::CheckSum($function);
        my $role = $checksum{$check};
        print "Copying $subDir/$name for role $role.\n";
        File::Copy::Recursive::fcopy("$dir/$subDir/$name.smp", "Profiles/$role.smp") || die "Copy failed for $role: $!";
    }
}
print "All done.\n";
