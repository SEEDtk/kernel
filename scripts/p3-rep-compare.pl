=head1 Display Statistics on Multiple Sets of Representative Genomes

    p3-rep-stats.pl [options] repDir1 repDir2 ... repDirN

This script produces a simple report comparing basic statistics for multiple representative-genome sets.
In particular, it displays the size of the largest group, the number of singletons, the similarity threshold,
the mean group size (excluding the largest and the singletons), the total number of groups, and the number
of genomes in groups with 100 or more members.


=head2 Parameters

Thepositional parameters are the names of the directories containing the databases. The files in the directories are
described in L<RepGenomeDb>.

The command-line options are as follows.

=over 4

=item verbose

If specified, status messages will be written to the standard error output.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use Math::Round;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir1 repDir2 ... repDirN',
        ['verbose|debug|v', 'display status messages on STDERR']);
# Get the input directories.
my (@repDirs) = @ARGV;
if (! @repDirs) {
    die "No input directories specified.";
} else {
    my @bad = grep { ! -d $_ } @repDirs;
    if (@bad) {
        die "Missing or invalid input directories: " . join(', ', @bad);
    }
}
# Get the options.
my $debug = $opt->verbose;
# These hashes contain the results.  All are keyed by directory name.
my (%maxGroup, %singles, %meanSize, %totGroups, %sim, %genomes, %hunGroups);
# Loop through the databases.
for my $repDir (@repDirs) {
    # Load the database.
    print STDERR "Loading database from $repDir.\n" if $debug;
    my $repDB = RepGenomeDb->new_from_dir($repDir);
    # Get the list of genomes.
    print STDERR "Counting set sizes.\n" if $debug;
    my $genomeList = $repDB->rep_list();
    # Get the similarity threshold.
    $sim{$repDir} = $repDB->score();
    # Compute the totals for this directory.
    my $hunGroups = 0;
    my $maxGroup = 1;
    my $singles = 0;
    my $totGroups = 0;
    my $totSize = 0;
    my $sizeGroups = 0;
    my $genomes = 0;
    for my $genomeID (@$genomeList) {
        my $repGenome = $repDB->rep_object($genomeID);
        my $repList = $repGenome->rep_list();
        my $count = scalar @$repList + 1;
        $totGroups++;
        $genomes += $count;
        if ($count == 1) {
            $singles++;
        } else {
            $sizeGroups++;
            $totSize += $count;
            if ($count > $maxGroup) {
                $maxGroup = $count;
            }
            if ($count >= 100) {
                $hunGroups += $count;
            }
        }
    }
    print STDERR "$singles singleton groups found out of $totGroups covering $genomes genomes.\n" if $debug;
    # Compute the interior mean.
    my $meanSize = 1;
    if ($maxGroup > 1) {
        $totSize -= $maxGroup;
        $sizeGroups--;
        if ($sizeGroups > 0) {
            $meanSize = Math::Round::nearest(0.01, $totSize / $sizeGroups);
        }
    }
    $maxGroup{$repDir} = $maxGroup;
    $meanSize{$repDir} = $meanSize;
    $singles{$repDir} = $singles;
    $totGroups{$repDir} = $totGroups;
    $genomes{$repDir} = $genomes;
    $hunGroups{$repDir} = $hunGroups;
}
# Now produce the output.
print STDERR "Writing output.\n" if $debug;
P3Utils::print_cols([qw(name similarity genomes groups singles ge100 largest mean_size)]);
for my $repDir (@repDirs) {
    P3Utils::print_cols([$repDir, $sim{$repDir}, $genomes{$repDir}, $totGroups{$repDir}, $singles{$repDir},
            $hunGroups{$repDir}, $maxGroup{$repDir}, $meanSize{$repDir}]);
}
print STDERR "All done.\n" if $debug;