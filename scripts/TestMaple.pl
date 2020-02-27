use strict;
use FIG_Config;

open(my $ih, '<', "subRoles.tbl") || die "Could not open subRoles.tbl: $!";
my %roles;
my $line = <$ih>;
while (! eof $ih) {
    if (<$ih> =~ /\t(.+)$/) {
        $roles{$1} = 1;
    }
}
close $ih; undef $ih;
my ($count, $kept) = 0;
my %fams;
open($ih, '<', '/vol/patric3/fams/2018-0531/merge1/prop.out/merged.families.1.1.nr-only') || die "Could not open families file: $!";
while (! eof $ih) {
    my $line = <$ih>;
    $count++;
    my ($fam, undef, undef, undef, undef, $role) = split /\t/, $line;
    if ($roles{$role}) {
        $fams{$fam} = 1;
        $kept++;
    }
    if ($count % 5000 == 0) {
        print STDERR "$count records process, $kept counted.\n";
    }
}
print scalar(keys %fams) . " families found for " . scalar(keys %roles) . " subsystem roles in $kept records.\n";
