#!/usr/bin/env perl
=head1 GTO Role Extraction

    gto_to_roles.pl [options] gtoFile roles.in.subsystems

This program extracts the roles from a L<GenomeTypeObject> and computes the role IDs.
It produces a list of warning messages about the unmapped roles on the standard error
output. The standard output itself is a three-column tab-delimited file. The first column
is a role name, the second column is a role ID, and the third column is a feature ID.

=head2 Parameters

The positional parameters are the name of the GTO input file (in json format) and the
name of a file containing the IDs of roles to keep (produced by L<build_role_tables.pl>).

The following command-line options are supported.

=over 4

=item counts

If specified, then the standard output will be a two-column file, each record consisting of
of a role ID followed by a count of the number of times the role was found in the input file.
In this case, the normal standard output is suppressed.

=cut

use strict;
use warnings;
use Data::Dumper;

use GenomeTypeObject;
use SeedUtils;
use ScriptUtils;
use RoleParse;

use Data::Dump qw(pp);
# die pp(\%INC);
my $opt = ScriptUtils::Opts('gto roles.in.subsystems',
        ['counts|count|k', 'produce counts table']);
my ($gto, $roles_in_subsystems) = @ARGV;
die "Input file '$gto' does not exist" unless (-s $gto);
die "Input file '$roles_in_subsystems' does not exist" unless (-s $roles_in_subsystems);
my $count = $opt->counts;
my %roleMap = map { my ($roleID, $checksum) = split /\t/;
                         ($checksum => $roleID)
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
        my $checksum = RoleParse::Checksum($role);
        # Compute the ID for this checksum.
        my $id = $roleMap{$checksum};
        if ($id) {
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
