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

=head1 Build Role Tables: Function Predictors Step 1

    build_role_tables.pl [ options ] outDir

This script builds two tab delimited files in a specified output directory. It is the first step in creating
function predictors. The second step is L<build_matrix.pl>.

=over 4

=item raw.table

This file contains (0) a role description, (1) a role ID, and (2) a feature ID in each record. There is one record for
every feature in well-behaved genomes that is currently in a subsystem.

=item roles.in.subsystems

This file contains (0) a role ID, (1) a role checksum, and (2) a role description in each record. There is one record for every
role appearing in the above file.

=back

=head2 Parameters

The positional parameter is the name of the directory to contain the output files.

The command-line options are those found in L<Shrub/script_options>.
=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
                Shrub::script_options(),
        );
# Create a statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Verify the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Invalid output directory $outDir.";
}
# This hash prevents duplicates. It maps each role ID found to a sub-hash keyed by feature ID.
my %rolePegs;
# This hash maps each role ID to its description. We populate it with the roles in core subsystems.
print "Computing subsystem roles.\n";
my %roles = map { $_->[0] => [$_->[1], $_->[2]] } $shrub->GetAll('Subsystem Role', 'Subsystem(privileged) = ?', [1], 'Role(id) Role(checksum) Role(description)');
$stats->Add(subsysRoles => scalar keys %roles);
# Open the main output file.
open(my $oh, ">$outDir/raw.table") || die "Could not open raw.table: $!";
# Loop through the desired features.
print "Looking for useful pegs.\n";
my $q = $shrub->Get('Role2Function Function2Feature Feature Feature2Genome Genome',
        'Function2Feature(security) = ? AND Feature(feature-type) = ? AND Genome(well-behaved) = ?', [2, 'peg', 1],
        'Role2Function(from-link) Function2Feature(to-link)');
print "Looping through useful pegs.\n";
while (my $roleData = $q->Fetch()) {
    my ($roleID, $peg) = $roleData->Values(['Role2Function(from-link)', 'Function2Feature(to-link)']);
    my $processed = $stats->Add(pegsFound => 1);
    # Only proceed if this is a role in a subsystem.
    if (! $roles{$roleID}) {
        $stats->Add(pegNoSubsysRole => 1);
    } elsif ($rolePegs{$roleID}{$peg}) {
        $stats->Add(pegDuplicate => 1);
    } else {
        # Output this role/peg pair.
        print $oh "$roles{$roleID}[1]\t$roleID\t$peg\n";
        $rolePegs{$roleID}{$peg} = 1;
        $stats->Add(pegOutput => 1);
    }
    print "$processed pegs processed.\n" if ($processed % 10000 == 0);
}
# Now output the list of roles used.
close $oh; undef $oh;
open($oh, ">$outDir/roles.in.subsystems") || die "Could not open roles.in.subsystems: $!";
for my $roleID (sort keys %rolePegs) {
    print $oh "$roleID\t$roles{$roleID}[0]\t$roles{$roleID}[1]\n";
    $stats->Add(roleOutput => 1);
}
print "All done.\n" . $stats->Show();
