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
use RepoLoader;
use RoleParse;
use SeedUtils;
use Stats;

=head1 Produce a Report of Role IDs

    analyze_role_ids.pl [ options ] inputDir

This script runs through all the assignments in the genome repository and checks for role synonyms, that is, roles with different
strings but the same ID. It also lists roles that have similar IDs, in case they should be synonoms.

=head2 Parameters

The positional parameter is the name of the L<genome repository directory|ExchangeFormat>.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inputDir',
                Shrub::script_options(),
        );
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the input repository.
my ($inputDir) = @ARGV;
if (! $inputDir) {
    die "No repository directory specified.";
} elsif (! -d $inputDir) {
    die "Invalid directory name $inputDir.";
}
my $loader = RepoLoader->new();
# Get the hash of genomes to directory names.
my $gHash = $loader->FindGenomeList($inputDir);
# This is a two-level hash that will map role IDs to role descriptions and then to number of feature IDs.
my %roleMap;
# This hash will map role checksums to role IDs.
my %roleIDs;
# This hash will gather similar IDs.
my %roleGroups;
# Create a hash of role IDs to subsystem names.
my %roleSS;
my $q = $shrub->Get('Role2Subsystem Subsystem', '', [], 'Role2Subsystem(from-link) Subsystem(name)');
while (my $ssDatum = $q->Fetch()) {
    my ($roleID, $ssName) = $ssDatum->Values(['Role2Subsystem(from-link)', 'Subsystem(name)']);
    $roleSS{$roleID}{$ssName} = 1;
}
for my $roleID (keys %roleSS) {
    my $subHash = $roleSS{$roleID};
    $roleSS{$roleID} = [sort keys %$subHash];
}
# Loop through the genomes.
for my $genome (sort keys %$gHash) {
    my $gData = $gHash->{$genome};
    my ($dir, $name) = @$gData;
    $stats->Add(genomeIn => 1);
    print STDERR "Processing $genome: $name.\n";
    my @fidFiles = map { "$dir/$_" } qw(non-peg-info peg-info);
    for my $fidFile (@fidFiles) {
        $stats->Add(fileIn => 1);
        open(my $ih, '<', $fidFile) || die "Could not open feature file $fidFile: $!";
        while (! eof $ih) {
            my $line = <$ih>;
            $stats->Add(lineIn => 1);
            my ($fid, $function) = ($line =~ /^(fig\|\d+\.\d+\.\w+\.\d+)\s\S+\s(.+)/);
            $stats->Add(featureIn => 1);
            if (! $fid) {
                $stats->Add(badLine => 1);
            } else {
                my @roles = SeedUtils::roles_of_function($function);
                $stats->Add(functionIn => 1);
                for my $role (@roles) {
                    my $checkSum = RoleParse::Checksum($role);
                    $stats->Add(roleIn => 1);
                    if (! exists $roleIDs{$checkSum}) {
                        # Here we are seeing the role for the first time.
                        my ($roleID) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checkSum], 'id');
                        if (! $roleID) {
                            # The role does not exist in Shrub. Blank out the checksum's record.
                            $stats->Add(roleNotFound => 1);
                            $roleIDs{$checkSum} = '';
                        } else {
                            # The role exists. Save the ID.
                            $roleIDs{$checkSum} = $roleID;
                            # Separate out the prefix and suffix to put the role in a group.
                            if ($roleID =~ /^(.+?)\d+$/) {
                                $roleGroups{$1}{$roleID} = 1;
                                $stats->Add(suffixedRole => 1);
                            } else {
                                $roleGroups{$roleID}{$roleID} = 1;
                                $stats->Add(unsuffixedRole => 1);
                            }
                        }
                    } else {
                        # Here we have seen the role before.
                        $stats->Add(oldRoleID => 1);
                    }
                    # Now we are sure to have something in the role ID hash.
                    my $roleID = $roleIDs{$checkSum};
                    if ($roleID) {
                        # It's a real role. Save this text and feature ID.
                        $roleMap{$roleID}{$role}++;
                        $stats->Add(fidRoleStored => 1);
                    }
                }
            }
        }
    }
}
# Now all the genomes have been processed. Look for role IDs with multiple descriptions.
for my $roleID (sort keys %roleMap) {
    my $subHash = $roleMap{$roleID};
    my @descs = sort keys %$subHash;
    if (scalar(@descs) <= 1) {
        $stats->Add(roleSingleton => 1);
    } else {
        $stats->Add(roleWithSynonyms => 1);
        print "$roleID\n";
        my $subs = $roleSS{$roleID} // [];
        for my $sub (@$subs) {
            print "\t** IN:\t$sub\n";
        }
        for my $desc (@descs) {
            print "\t$desc\t$subHash->{$desc}\n";
        }
    }
}
print "\n\n";
# Now look for groups of similar role IDs.
for my $roleGroup (sort keys %roleGroups) {
    my @subRoles = sort keys %{$roleGroups{$roleGroup}};
    if (scalar(@subRoles) > 1) {
        $stats->Add(multiRoleGroup => 1);
        print "$roleGroup\n";
        for my $subRole (@subRoles) {
            my ($desc) = sort keys %{$roleMap{$subRole}};
            print "\t$subRole\t$desc\n";
        }
    }
}
print "\n\n";
print STDERR "All done.\n" . $stats->Show();