#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use GenomeTypeObject;
use SeedUtils;
use ScriptUtils;
use Shrub;
use Shrub::Roles;

use Data::Dump qw(pp);
# die pp(\%INC);
my $opt = ScriptUtils::Opts('gto roles.in.subsystems', Shrub::script_options(),
        ['counts|count|k', 'produce counts table']);
my ($gto, $roles_in_subsystems) = @ARGV;
my $shrub = Shrub->new_for_script();
die "Input file '$gto' does not exist" unless (-s $gto);
die "Input file '$roles_in_subsystems' does not exist" unless (-s $roles_in_subsystems);
my $count = $opt->counts;
my %safeRoles = map { my ($roleID) = split /\t/;
                         ($roleID => 1)
} &SeedUtils::file_read($roles_in_subsystems);
# die Dumper(\%roles_to_IDs);


my $proc = GenomeTypeObject->new({file => $gto});
# die Dumper($proc);

my @CDSs = map {
    [ $_->{id}, $_->{function} ]
} grep {
    ($_->{type} eq q(CDS)) || ($_->{id} =~ m{\.peg\.})
} $proc->features();
my %counts;
foreach my $cds (@CDSs) {
    my ($fid, $func) = @$cds;
    my @roles = &SeedUtils::roles_of_function($func);
    foreach my $role (@roles) {
        # Compute the role's checksum.
        my $checksum = Shrub::Roles::Checksum($role);
        # Compute the ID for this checksum.
        my ($id) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checksum], 'id');
        if ($id && $safeRoles{$id}) {
            if ($count) {
                $counts{$id}++;
            } else {
                print STDOUT (join("\t", ($role, $id, $fid)), "\n");
            }
        }
        else {
            warn "Could not map role:\t$fid\t$role\n";
        }
    }
}
if ($count) {
    for my $roleID (sort keys %counts) {
        print join("\t", $roleID, $counts{$roleID}) . "\n";
    }
}
