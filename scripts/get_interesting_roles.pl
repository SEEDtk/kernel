#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


use strict;
use warnings;
use FIG_Config;
use Shrub;
use ScriptUtils;
use Stats;

=head1 Create File of Roles in Subsystems

    get_interesting_roles.pl [ options ] fileName

This script creates a simple tab-delimited file of the roles found in subsystems (that is, the roles that are best defined,
and therefore most interesting). This file can then be used for role processing in P3 scripts.

The output file will have one record per role, and each record will contain (0) the role ID, (1) the role checksum, and (2) the
role description.

=head2 Parameters

The single positional parameter is the name for the output file.

The command-line options are those found in L<Shrub/script_options>.

=cut

$| = 1;
# Create the statistics object.
my $stats = Stats->new();
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('fileName', Shrub::script_options());
# Get the output file name.
my ($fileName) = @ARGV;
if (! $fileName) {
    die "No output file name specified.";
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# This hash will contain the output. It is keyed by role ID, and maps each role to a [checksum,description] pair.
my %roles;
# Pull the roles out of the database. We do the hash thing because if a role is in two subsystems we'll see it twice.
print "Reading roles from database.\n";
my $qh = $shrub->Get("Role Role2Subsystem", '', [], [qw(id checksum description)]);
while (my $resultRow = $qh->Fetch()) {
    $stats->Add(roleRow => 1);
    my $id = $resultRow->PrimaryValue('id');
    my $checksum = $resultRow->PrimaryValue('checksum');
    my $description = $resultRow->PrimaryValue('description');
    $roles{$id} = [$checksum, $description];
}
# Open the output file.
print "Writing roles to output.\n";
open(my $oh, '>', $fileName) || die "Could not open output file $fileName: $!";
# Write out the roles found.
for my $id (sort keys %roles) {
    $stats->Add(roleOut => 1);
    print $oh join("\t", $id, @{$roles{$id}}) . "\n";
}
close $oh;
print "All done.\n" . $stats->Show();
