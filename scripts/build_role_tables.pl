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

This script builds three tab delimited files in a specified output directory. It is the first step in creating
function predictors. The second step is L<build_matrix.pl>.

IMPORTANT NOTE: /vol/patric3/fams/I<date>/kmers/function.index contains the functions currently in PATRIC.

=over 4

=item raw.table

This file contains (0) a role description, (1) a role ID, and (2) a feature ID in each record. There is one record for
every feature in well-behaved genomes that is currently in a subsystem.

=item roles.in.subsystems

This file contains (0) a role ID, (1) a role checksum, and (2) a role description in each record. There is one record for every
role appearing in the above file.

=item roles.to.use

This file contains (0) a role ID. There is one record for every role that appears in at least 100 genomes.

=back

=head2 Parameters

The positional parameter is the name of the directory to contain the output files.

The command-line options are those found in L<Shrub/script_options>.

=head2 Basic Procedure for Creating the Predictors.

First, you must identify a directory for the predictors and a directory for the role files. In the example below,
we will use C<FunctionPredictors.X> for the predictors and C<Staging> for the role files.

The first command creates the initial role files.

    build_role_tables Staging

The following four commands are run repeatedly in a loop. When the output of L<build_roles_to_use.pl> indicates that
all of the roles are good, the loop stops.

    build_matrix --clear Staging/raw.table FunctionPredictors.X Staging/roles.to.use
    build_LDA_predictors FunctionPredictors.X
    cp Staging/roles.to.use FunctionPredictors.X/
    build_roles_to_use FunctionPredictors.X Staging

Once the role set has converged, to install the new predictors you must (1) copy C<roles.in.subsystems> from C<Staging>
to the SEEDtk Global directory. (2) Copy C<FunctionPredictors.X/roles.to.use> to the SEEDtk Global directory. (3)
Symbolically link C<FunctionPredictors> in the global directory to C<FunctionPredictors.X>.

=cut

$| = 1;
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
print "Writing roles.in.subsystems and roles.to.use.\n";
open($oh, ">$outDir/roles.in.subsystems") || die "Could not open roles.in.subsystems: $!";
open(my $rh, ">$outDir/roles.to.use") || die "Could not open roles.to.use: $!";
for my $roleID (sort keys %rolePegs) {
    print $oh "$roleID\t$roles{$roleID}[0]\t$roles{$roleID}[1]\n";
    $stats->Add(roleOutput => 1);
    my @pegs = keys %{$rolePegs{$roleID}};
    my %genomes;
    for my $peg (@pegs) {
        my $genome = SeedUtils::genome_of($peg);
        $genomes{$genome} = 1;
    }
    my $count = scalar keys %genomes;
    if ($count >= 100) {
        print "$count genomes found for $roleID.\n";
        print $rh "$roleID\n";
        $stats->Add(roleUsed => 1);
    }
}
print "All done.\n" . $stats->Show();
