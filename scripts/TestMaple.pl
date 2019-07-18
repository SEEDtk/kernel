use strict;
use FIG_Config;
use SeedUtils;
use P3Utils;

opendir(my $dh, 'Bins_HMP') || die "Could not open binning directory: $!";
my @samples = grep { -s "Bins_HMP/$_/Eval/index.tbl" } readdir $dh;
closedir $dh;
my $binCount = 0;
my %roleBad;
for my $sample (@samples) {
    print STDERR "Processing $sample.\n";
    my $subDir = "Bins_HMP/$sample";
    # Loop through the bins.
    for (my $i = 1; -s "$subDir/bin$i.gto"; $i++) {
        my $gto = SeedUtils::read_encoded_object("$subDir/bin$i.gto");
        my $qHash = $gto->{quality};
        # We only process bins that Maulik thinks are good.
        if ($qHash->{genome_quality} eq 'Good') {
            print STDERR "Scanning $sample bin $i.\n";
            $binCount++;
            my %roles;
            for my $type (qw(consistency_roles completeness_roles)) {
                my $rCounts = $qHash->{problematic_roles_report}{$type};
                for my $role (keys %$rCounts) {
                    if ($rCounts->{$role}[0] != $rCounts->{$role}[1]) {
                        $roles{$role} = 1;
                    }
                }
            }
            for my $role (keys %roles) {
                $roleBad{$role}++;
            }
        }
    }
}
print STDERR "$binCount bins scanned.\n";
if ($binCount > 0) {
    my @roles = sort { $roleBad{$b} <=> $roleBad{$a} } keys %roleBad;
    print "Role\tcount\tpercent\n";
    for my $role (@roles) {
        my $pct = $roleBad{$role} * 100 / $binCount;
        P3Utils::print_cols([$role, $roleBad{$role}, $pct]);
    }
}