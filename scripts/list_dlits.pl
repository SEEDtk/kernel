#!/usr/bin/env perl
use FIG_Config;
my $wrong_env = 1;
eval {
    require FIG;
    $wrong_env = 0;

};
if ($wrong_env) {
    die "Invalid environment.  THis must be run in CoreSEED mode.";
}
my $fig = FIG->new();
print STDERR "Requesting DLITs.\n";
my $dlits = $fig->all_dlits();
print STDERR scalar(@$dlits) . " DLITs found.\n";
print "status\taa_sequence_md5\tpubmed\tcurator\tGO_code\n";
for my $dlit (@$dlits) {
    print join("\t", @$dlit) . "\n";
}
