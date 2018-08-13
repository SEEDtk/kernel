=head1 Analyze Universal Roles By Taxonomy

    p3-taxon-analysis.pl [options] workDir

This script determines the universal roles for taxonomic groupings. A I<universal role> is one that occurs singly 97% of the time in genomes
of a given taxonomic grouping. The script takes as input a file of PATRIC genomes believed to complete along with their singly-occurring roles
and taxonomic lineages.

Progress messages are sent to the standard output. The work directory must contain the file C<genomes.tbl> produced by L<p3-genome-unis.pl>. The
file C<roles.tbl> will be created to contain a list of the useful taxonomic groupings, one per line, with the taxonomic ID, the name, the
number of genomes in the group, and a list of the required role IDs. This can be fed into L<group_marker_roles.pl> and L<compute_taxon_map.pl>
to finish building the EvalG directory structure.

=back

=head2 Parameters

The positional parameter is the name of the work directory. The work directory must contain the C<genomes.tbl> file produced by
L<p3-genome-unis.pl>.

The options in L<Shrub/script_options> can be used to specify the Shrub database. This database contains the full taxonomy tree.

Additional options include the following.

=over 4

=item min

Minimum percentage of genomes in a taxonomic grouping that must contain a role for it to be considered universal. The default is C<97>.

=item size

Minimum number of genomes in a taxonomic grouping required for the grouping to be considered worthwhile. The default is C<100>.

=item rMin

Minimum number of roles for a taxonomic grouping for it to be considered useful in determining completeness. The default is C<60>.

=item strict

If specified, the universal roles for a taxonomic grouping will be computed by taking the intersection of the universal roles for the lower groupings.
Otherwise, the lower groupings will be weighted by size.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use GEO;
use File::Copy::Recursive;
use EvalCon;
use Stats;
use Shrub;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('workDir', Shrub::script_options(),
        ['min|m=f', 'minimum threshold percentage', { default => 97 }],
        ['size|s=i', 'minimum number of genomes per group', { default => 100 }],
        ['rMin|r=i', 'minimum number of roles for a useful group', { default => 60 }],
        ['strict', 'use strict merging']
        );
my $stats = Stats->new();
# Get the output directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No output directory specified.";
} elsif (! -d $workDir) {
    print "Creating $workDir.\n";
    File::Copy::Recursive::pathmk($workDir) || die "Could not create $workDir: $!";
}
print "Accessing data.\n";
# Get access to Shrub.
my $shrub = Shrub->new_for_script($opt);
# Open the input file and skip the header.
open(my $ih, "<$workDir/genomes.tbl") || die "Could not open genomes.tbl: $!";
my $line = <$ih>;
# This hash will contain the number of genomes in each taxonomic group found.
my %taxSize;
# This hash will map each taxonomic grouping to a sub-hash of the number of times each role occurs singly.
my %taxRoles;
# This hash maps each taxonomic grouping to its parent.
my %parent;
# This hash tracks the groups that have sub-groups.
my %parental;
# Here we get the size and count of each taxonomic group. We also memorize the parent-child relationships.
my $count = 0;
my $start = time;
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($genome, $lineageX, $rolesX) = split /\t/, $line;
    my @lineage = split /,/, $lineageX;
    my @unis = split /,/, $rolesX;
    # Compute the universal roles.
    $stats->Add(genomesProcessed => 1);
    $stats->Add(unisFound => scalar @unis);
    # Now sum us into all the taxonomic groups!
    for my $taxon (@lineage) {
        $stats->Add(genomeSummed => 1);
        for my $role (@unis) {
            $taxRoles{$taxon}{$role}++;
        }
        $taxSize{$taxon}++;
    }
    # Create the lineage connections.
    my $child = pop @lineage;
    while (my $parent = pop @lineage) {
        # Record this parent-child relationship.
        $parent{$child} = $parent;
        $parental{$parent} = 1;
        $child = $parent;
    }
    $count++;
    if ($count % 1000 == 0) {
        my $duration = time - $start;
        print STDERR "$count genomes processed in $duration seconds.\n";
    }
}
# We are ready to look for groups to output. Open the output file.
open(my $oh, ">$workDir/roles.tbl") || die "Could not open roles.tbl: $!";
P3Utils::print_cols(['Taxon', 'Name', 'Size', 'Roles'], oh => $oh);
my %outGroup;
# Now we process the leaf groups. If we are strict, then roles that do not qualify are cleared from the parent groups all the way
# up the chain.
my $strict = $opt->strict;
print "Processing leaf groups in " . ($strict ? 'strict' : 'statistical') . " mode.\n";
my $minF = $opt->min;
my $rMin = $opt->rmin;
$count = 0;
for my $taxon (grep { ! $parental{$_} } sort keys %taxRoles) {
    $stats->Add(leafGroup => 1);
    ProcessTaxon($taxon, $minF, $strict);
    $count++;
    print STDERR "$count leaf groups processed.\n" if $count % 1000 == 0;
}
# Now we process the other groups. Note that if we are in strict mode, roles that were too infrequent in lower groups will have a
# count of 0 and will automatically fail.
print "Processing parent groups.\n";
$count = 0;
for my $taxon (sort keys %parental) {
    $stats->Add(parentGroup => 1);
    ProcessTaxon($taxon, $minF, 0);
    $count++;
    print STDERR "$count parent groups processed.\n" if $count % 1000 == 0;
}
# All done.
print "All done.\n" . $stats->Show();

## Find the acceptable roles in this taxonomic grouping. Output the group if it qualifies.
sub ProcessTaxon {
    my ($taxon, $minF, $strict) = @_;
    # Compute the threshold for this group. Note the trick to emulate ceil(x), which works because the denominator is 100.
    my $limit = int($taxSize{$taxon} * $minF / 100 + 0.99);
    # Find roles that exceed the threshold.
    my $roleCounts = $taxRoles{$taxon};
    # This will contain the roles found.
    my @found;
    for my $role (sort keys %$roleCounts) {
        if ($roleCounts->{$role} >= $limit) {
            push @found, $role;
            $stats->Add(roleCandidate => 1);
        }
    }
    # Get the size of the group, the role count, and the name.
    my $groupSize = $taxSize{$taxon};
    my $groupRoles = scalar(@found);
    my ($taxName) = $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(id) = ?', [$taxon], 'scientific-name');
    if (! $taxName) {
        # A missing name is serious. It means the group is messed up.
        $stats->Add(nameNotFound => 1);
        $taxName = 'unknown';
    }
    if ($groupRoles < $rMin) {
        $stats->Add(groupTooFewRoles => 1);
        print "Group $taxon ($taxName) has only $groupRoles eligible roles.\n";
    } else {
        # Too few roles means we ignore the group. For other cases, we still propagate the group's statistics.
        if ($groupSize < $opt->size) {
            $stats->Add(groupTooSmall => 1);
        } elsif ($taxName) {
            # This group is worthwhile. It's big enough, has enough roles, and we found the name.
            P3Utils::print_cols([$taxon, $taxName, $groupSize, @found], oh => $oh);
            print "$taxon ($taxName) has $groupSize genomes and $groupRoles universal roles. Saved.\n";
            $stats->Add(groupStored => 1);
        }
        # If we are strict, propagate the group's rejected roles.
        if ($strict) {
            # Create a hash of the eligible roles in this group.
            my %goodRoles = map { $_ => 1 } @found;
            # Loop up through the parents.
            my $pTaxon = $parent{$taxon};
            while ($pTaxon) {
                $stats->Add(strictGroupCheck => 1);
                # Zero out role counts for ineligible roles.
                my $pRoleCounts = $taxRoles{$pTaxon};
                for my $role (keys %$pRoleCounts) {
                    if (! $goodRoles{$role}) {
                        $pRoleCounts->{$role} = 0;
                        $stats->Add(strictRoleDeleted => 1);
                    }
                }
                $pTaxon = $parent{$pTaxon};
            }
        }
    }
}
