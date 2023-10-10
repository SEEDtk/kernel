#!/usr/bin/env perl
use FIG_Config;
my $wrong_env = 1;
eval {
    require FIG;
    $wrong_env = 0;

};
if ($wrong_env) {
    die "Invalid environment.  This must be run in CoreSEED mode.";
}
my $fig = FIG->new();
print STDERR "Requesting DLITs.\n";
my $dlits = $fig->all_dlits();
print STDERR scalar(@$dlits) . " DLITs found.\n";
print "status\taa_sequence_md5\tpubmed\tcurator\tGO_code\n";
my $out = 0;
for my $dlit (@$dlits) {
    if ($dlit->[1] =~ /^[a-f0-9]+$/) {
        print join("\t", @$dlit) . "\n";
        $out++;
    }
}
print STDERR "$out DLITs printed.\n";
