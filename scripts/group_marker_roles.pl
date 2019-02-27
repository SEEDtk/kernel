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
use POSIX qw(ceil);
use Time::HiRes;
use Data::Dump;
use P3Utils;

=head1 Group Marker Roles and Assign Weights

    group_marker_roles.pl [ options ] workDir

This script reads the C<roles.tbl> file produced by L<p3-taxon-analysis.pl> or L<p3-rep-common-roles.pl> and organizes the roles
into cluster groups. Each pair of incoming roles is examined against the L<Shrub> database to see which are commonly clustered. For
each grouping in C<roles.tbl>, the roles are organized into groups based on the pairs, and a new C<weighted.tbl> file is
written with weights assigned to the roles. The output file will be tab-delimited with two line types. Each header line
will contain all of the data columns from the group line. Each detail line will contain a role ID and the weight of the role.
A double-slash (C<//>) will be used as a delimiter.

=head2 Parameters

The positional parameter is the name of the work directory containing the input C<roles.tbl>, and into which the C<weighted.tbl>
will be written.  The role list must begin in the column named C<Roles>.  The first two columns must be the group ID and name.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item min

The minimum percent co-occurrence required for a role to be considered clustered. The default is C<95>, indicating the
roles must appear together 95% of the time.

=back

=head2 Basic Algorithm

For each group, we sort the roles and then get the features for each role. For each feature, we count one for
each clustered feature's roles. A role that appears enough times is considered clustered. Once we know two roles are clustered,
we save the information so we don't have to check again.

When we have the complete list of pairs for the group, we use the transitive law to fully cluster them and produce the output.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir',
        Shrub::script_options(),
        ['min|m=f', 'minimum percentage for a valid colocation', { default => 95 }],
        );
my $stats = Stats->new();
# Get the parameters.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No work directory specified.";
} elsif (! -s "$workDir/roles.tbl") {
    die "$workDir is not a directory or does not contain a roles.tbl file.";
}
my $min = $opt->min;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# This will contain the known pairs. Each pair is identified by the two role IDs, sorted and space-delimited.
# A value of C<1> means a cluster, a value of C<-1> means not a cluster.
my %pairs;
# This will contain the count of features for each role processed.
my %roleFeats;
# Open the output file.
open(my $oh, ">$workDir/weighted.tbl") || die "Could not open weighted.tbl: $!";
# Open the input file and loop through the groups.
print "Work directory is $workDir.\n";
open(my $ih, "<$workDir/roles.tbl") || die "Could not open roles.tbl: $!";
my (undef, $cols) = P3Utils::find_headers($ih, roleFile => 'Roles');
my ($roleCol) = @$cols;
while (! eof $ih) {
    my $start = time;
    # Get the group's id, name, middle fields, and roles.
    my @fields = ScriptUtils::get_line($ih);
    my ($groupID, $name, @roles) = @fields;
    my @middle;
    for (my $shifted = 2; $shifted < $roleCol; $shifted++) {
        push @middle, shift @roles;
    }
    my $markerCount = scalar(@roles);
    print "Processing $groupID: $name. $markerCount marker roles.\n";
    # Sort the roles and save a copy. This is a two-pass process.
    my @sorted = sort @roles;
    @roles = @sorted;
    # Loop through the roles.
    while (my $role = shift @sorted) {
        # Eliminate all the role pairs we already know.
        my %residual;
        for my $role2 (@sorted) {
            if (! $pairs{"$role/$role2"}) {
                $residual{$role2} = 0;
            }
        }
        # Only proceed if we have unknown roles in the mix.
        if (! keys %residual) {
            $stats->Add(roleAllKnown => 1);
        } else {
            $stats->Add(roleChecked => 1);
            print ".";
            # Get the features for this role.
            my $fidCount = $roleFeats{$role};
            if (! $fidCount) {
                $fidCount = $shrub->GetCount('Role2Function Function2Feature', 'Role2Function(from-link) = ? AND Function2Feature(security) = ?',
                        [$role, 2], 'Function2Feature(to-link)');
                $roleFeats{$role} = $fidCount;
                $stats->Add(roleFeaturesRead => 1);
            } else {
                $stats->Add(roleFeaturesFound => 1);
            }
            # Now we count the number of pairings in %residual. Get all the clustered roles for each feature.
            my $q = $shrub->Get('Role2Function Function2Feature Feature2Cluster Cluster2Feature Feature2Function Function2Role',
                'Role2Function(from-link) = ? AND Function2Feature(security) = ? AND Feature2Function(security) = ?',
                [$role, 2, 2], 'Function2Feature(to-link) Function2Role(to-link)');
            while (my $record = $q->Fetch()) {
                my ($role2) = $record->PrimaryValue('Function2Role(to-link)');
                $residual{$role2}++;
            }
            my $minCount = int(($fidCount * $min + 99) / 100);
            # Determine the pairing of these roles.
            for my $role2 (keys %residual) {
                if ($residual{$role2} >= $minCount) {
                    $stats->Add(pairClustered => 1);
                    $pairs{"$role/$role2"} = 1;
                } else {
                    $stats->Add(pairIndependent => 1);
                    $pairs{"$role/$role2"} = -1;
                }
            }
        }
    }
    print "\n";
    # Now the status is known for all the role pairs in this taxonomic grouping. Form them into partitions.
    my $partitionList = partition(\@roles, \%pairs);
    my $groupCount = scalar(@$partitionList);
    print "$groupCount role groups found for $markerCount roles.\n";
    $stats->Add(groupDelta => ($markerCount - $groupCount));
    # Finally, we print this taxonomic group.
    print $oh join("\t", $groupID, $name, @middle) . "\n";
    $stats->Add(taxonOut => 1);
    for my $partition (@$partitionList) {
        $stats->Add(groupOut => 1);
        my $w =  1 / scalar keys %$partition;
        for my $role (keys %$partition) {
            print $oh "$role\t$w\n";
            $stats->Add(roleOut => 1);
        }
    }
    print $oh "//\n";
    print int(time - $start) . " seconds to process group.\n";
}
print "All done.\n" . $stats->Show();

## This method creates partitions from known pairing. It accepts as input a sorted list of role IDs and a hash containing
## 1 for pairs and -1 for non-pairs. All possible pairs in the list will be represented. The roles are put into groups such
## that every role is in the same group as a role with which it is paired. The role list is destroyed.
sub partition {
    my ($list, $pairHash) = @_;
    # This will be the return, a list of hashes. Each hash is a group.
    my @retVal;
    # This hash maps a role to its group, the group identified as a position in the list.
    my %roles;
    # Loop through the roles.
    while (my $role = shift @$list) {
        # Get this role's group ID and hash.
        my $group = $roles{$role};
        my $gHash;
        if (! defined $group) {
            $group = scalar(@retVal);
            $gHash = { $role => $group };
            $roles{$role} = $group;
            push @retVal, $gHash;
            $stats->Add(groupCreated => 1);
        } else {
            $gHash = $retVal[$group];
            $stats->Add(groupFound => 1);
        }
        for my $role2 (@$list) {
            # Is this other role paired with our role?
            if ($pairHash->{"$role/$role2"} == 1) {
                # Yes. Get the other role's group.
                my $group2 = $roles{$role2};
                if (! defined $group2) {
                    # The other role is not in a group. Add it to this one.
                    $gHash->{$role2} = 1;
                    $roles{$role2} = $group;
                    $stats->Add(groupRoleAdd => 1);
                } elsif ($group == $group2) {
                    # The other role is already in the same group.
                    $stats->Add(groupTransitive => 1);
                } else {
                    # Here we must merge the groups.
                    my $group2H = $retVal[$group2];
                    for my $roleX (keys %$group2H) {
                        $roles{$roleX} = $group;
                        $gHash->{$roleX} = 1;
                        $stats->Add(groupMerged => 1);
                    }
                    $retVal[$group2] = undef;
                }
            }
        }
    }
    # Compress and return the list of hashes.
    my $retVal = [ grep { $_ } @retVal ];
    return $retVal;
}