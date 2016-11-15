#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use GenomeTypeObject;
use SeedUtils;
use ScriptUtils;

use Data::Dump qw(pp);
# die pp(\%INC);
my $opt = ScriptUtils::Opts('gto roles.in.subsystems',
        ['counts|count|k', 'produce counts table']);
my ($gto, $roles_in_subsystems) = @ARGV;
die "Input file '$gto' does not exist" unless (-s $gto);
die "Input file '$roles_in_subsystems' does not exist" unless (-s $roles_in_subsystems);
my $count = $opt->counts;
my %roles_to_IDs = map { chomp;
                         my ($roleID, $role) = split /\t/;
                         ($role => $roleID)
} &SeedUtils::file_read($roles_in_subsystems);
# die Dumper(\%roles_to_IDs);


my $proc = GenomeTypeObject->new({file => $gto});
# die Dumper($proc);

my @CDSs = map {
    [ $_->{id}, $_->{function} ]
} grep {
    ($_->{type} eq q(CDS)) || ($_->{id} =~ m{\.peg\.})
} $proc->features();

foreach my $cds (@CDSs) {
    my ($fid, $func) = @$cds;
    my @roles = &SeedUtils::roles_of_function($func);
    my %counts;
    foreach my $role (@roles) {
        if (my $roleID = $roles_to_IDs{$role}) {
            if ($count) {
                $counts{$roleID}++;
            } else {
                print STDOUT (join("\t", ($role, $roleID, $fid)), "\n");
            }
        }
        else {
            warn "Could not map role:\t$fid\t$role\n";
        }
    }
    if ($count) {
        for my $roleID (sort keys %counts) {
            print join("\t", $roleID, $counts{$roleID}) . "\n";
        }
    }
}
