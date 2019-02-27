=head1 Compute the Common Universal Roles for Represented Groups

    p3-rep-common-roles.pl [options] outDir repDir1 repDir2 ... repDirN

This program will create the database for an L<EvalCom/Rep> completeness/contamination evaluator.  All such systems
rely on a weighted sum of roles.  This script takes as input the output of L<p3-uni-roles.pl> and produces a C<roles.tbl>
that can be passed into L<group_marker_roles.pl>.

For each genome grouping, the output file will contain the representative genome ID (which serves as the group ID),
the representative genome name, and then the common roles.

Once everything is read in, we run through the groups computing common roles and writing the output file.

=head2 Parameters

The positional parameters are the name of the output directory and the names of the representative-genome directories
to process.

The input should be the file produced by L<p3-uni-roles.pl>.  This is a tab-delimited file with headers.  The
first two columns contain the genome ID and the seed protein sequence.  The remaining columns contain the marker roles.

=over 4

=item input

The name of the input file.  The default is C<uni_roles.tbl> in the output directory.

=item percent

The minimum percent of genomes that must have a role for it to be considered common.  The default is C<97>.

=item size

The minimum number of genomes in a group for it to be considered useful.  The default is C<100>.

=item minRoles

The minimum number of roles in a group for it to be considered useful.  The default is C<80>.

=back

=cut

use strict;
use P3Utils;
use RepGenomeDb;
use Stats;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir repDir1 repDir2 ... repDirN', P3Utils::ih_options(),
        ['percent|pct|p=i', 'percent threshold for common roles', { default => 97 }],
        ['size|sz|S=i', 'minimum size for a group to be interesting', { default => 100 }],
        ['input=s', 'name of input file (default uni_roles.tbl in output directory)'],
        ['minroles|minRoles|min=i', 'minimum number of common roles for a group to be interesting', { default => 80 }]
        );
my $stats = Stats->new();
# Get the options.
my $pct = $opt->percent;
my $size = $opt->size;
my $minRoles = $opt->minroles;
# Get the output directory and the representative-genome directories.
my ($outDir, @repDirs) = @ARGV;
# Check the representative-genome directory.
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
}
# Now we read the representative-genome directories.  For each genome, we want to know all its representatives.
# The hash below maps the directory's similarity score to the database object itself.
my %repDBs;
print "Reading representative-genome databases.\n";
for my $repDir (@repDirs) {
    my $repDB = RepGenomeDb->new_from_dir($repDir, verbose => 1, unconnected => 1);
    $repDBs{$repDB->score} = $repDB;
    $stats->Add(repDBIn => 1);
}
# Save the list of scores.
my @scores = sort keys %repDBs;
# Compute the input file.
print "Opening input file.\n";
my $inputFile = $opt->input // "$outDir/uni_roles.tbl";
# Open it and skip the header.
open(my $ih, '<', $inputFile) || die "Could not open $inputFile: $!";
my $line = <$ih>;
# Each representative genome corresponds to a possible group.  The group ID is the genome ID plus a slash and the directory's
# similarity score.  The group name is the genome name followed by the similarity score in square brackets with an R prefix
# (e.g. "[R50]").  These hashes are all keyed on group ID.  The first maps to the group name, the second to the number of
# genomes in the group, the third to a sub-hash that counts the number of times each role is found, and the fourth has the
# seed protein.
my (%groupName, %groupSize, %groupRoles, %groupSeed);
# Now we process the incoming genomes one at a time.
my $gCount = 0;
while (! eof $ih) {
    my $line = <$ih>;
    my ($genomeID, $prot, @roles) = P3Utils::get_fields($line);
    $stats->Add(genomeIn => 1);
    $stats->Add(rolesIn => scalar @roles);
    $gCount++;
    print "Processing $genomeID ($gCount).\n";
    # Use the protein sequence to find all the groups containing this genome.
    for my $score (@scores) {
        my $repDB = $repDBs{$score};
        my $repList = $repDB->list_reps($prot, $score);
        for my $rep (keys %$repList) {
            # Compute the group ID.
            my $groupID = init_group($rep, $score, $repDB, \%groupName, \%groupSeed);
            $groupSize{$groupID}++;
            $stats->Add(groupFound => 1);
            for my $role (@roles) {
                $groupRoles{$groupID}{$role}++;
                $stats->Add(roleRecorded => 1);
            }
        }
    }
}
# Now we have stored all of the information about all the genomes.  We loop through the groups, keeping the ones with
# adequate information.  First, we initialize the output file.
open(my $oh, '>', "$outDir/roles.tbl") || die "Could not open output file in $outDir: $!";
P3Utils::print_cols([qw(Group Name Prot Score Roles)], oh => $oh);
my $groupTot = scalar keys %groupSize;
print "Processing $groupTot groups.\n";
my $groupCount = 0;
for my $group (sort keys %groupSize) {
    $groupCount++;
    my $gSize = $groupSize{$group};
    if ($gSize < $size) {
        $stats->Add(groupTooSmall => 1);
    } else {
        # Here the group is reasonable.  Check its roles.
        my $minimum = $pct * $gSize / 100;
        print "Processing $group ($groupCount of $groupTot).\n";
        my $roleH = $groupRoles{$group};
        my @kept;
        for my $role (keys %$roleH) {
            if ($roleH->{$role} >= $minimum) {
                push @kept, $role;
                $stats->Add(roleKept => 1);
            } else {
                $stats->Add(roleRejected => 1);
            }
        }
        if (scalar @kept < $minRoles) {
            $stats->Add(groupFewRoles => 1);
        } else {
            # Here we are keeping this group.
            my (undef, $score) = split /\//, $group;
            $stats->Add(groupKept => 1);
            P3Utils::print_cols([$group, $groupName{$group}, $groupSeed{$group}, $score, sort @kept], oh => $oh);
        }
    }
}
print "All done.\n" . $stats->Show();

## Initialize a group:  insure it has a name and compute its ID.
sub init_group {
    my ($rep, $score, $repDB, $groupNameH, $groupSeedH) = @_;
    my $retVal = "$rep/$score";
    if (! exists $groupNameH->{$retVal}) {
        print "Creating group $retVal.\n";
        $stats->Add(groupCreated => 1);
        my $repO = $repDB->rep_object($rep);
        my $name = $repO->name();
        $groupNameH->{$retVal} = "$name [R$score]";
        $groupSeedH->{$retVal} = $repO->prot();
    }
    return $retVal;
}