=head1 Compute the Common Universal Roles for Represented Groups

    p3-rep-common-roles.pl [options] repDir

This program will analyze the represented-genome groups in a representative-genome directory.  For each represented set above a certain
size, it will compute the universal roles in common and write them to the output.

=head2 Parameters

The positional parameter is the name of the representative-genome directory.

Additional command-line options are those in L<EvalCon/role_options> plus the following.

=over 4

=item verbose

Display progress messages on STDERR.

=item resume

Resume after a previous failed run.  The parameter should be the ID of the last representative genome successfully processed.

=item percent

The minimum percent of genomes that must have a role for it to be considered common.  The default is C<97>.

=item size

The minimum number of genomes in a group for it to be considered useful.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use EvalCon;
use Stats;
use GEO;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('repDir', EvalCon::role_options(),
        ['verbose|debug|v', 'display progress messages on STDERR'],
        ['resume=s', 'resume after specified genome'],
        ['percent|pct|p=i', 'percent threshold for common roles', { default => 97 }],
        ['size|sz|S=i', 'minimum size for a group to be interesting', { default => 100 }],
        );
# Get the options.
my $debug = $opt->verbose;
my $logH = ($debug ? \*STDERR : undef);
my $resume = $opt->resume;
my $pct = $opt->percent;
my $size = $opt->size;
# Get access to PATRIC.
print STDERR "Connecting to PATRIC.\n" if $debug;
my $p3 = P3DataAPI->new();
# Check the representative-genome directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Input directory $repDir not found or invalid.";
}
# Read the role definitions.
my $evalCon = EvalCon->new_for_script($opt, $logH);
my ($nMap, $cMap) = $evalCon->roleHashes;
my $stats = $evalCon->stats;
# Compute the GEO options.
my %gOpts = (roleHashes => [$nMap, $cMap], logH => $logH, detail => 0, p3 => $p3, stats => $stats);
# Load it into memory.
print STDERR "Loading from $repDir.\n" if $debug;
my $repDB = RepGenomeDb->new_from_dir($repDir);
# Get the main list of representative genomes.
my $genomes = $repDB->rep_list();
my $totalGroups = scalar @$genomes;
my $nGroups = 0;
print STDERR "$totalGroups groups found in directory.\n" if $debug;
# Write the output file header.
if (! $opt->resume) {
    P3Utils::print_cols([qw(group name size minDist meanDist maxDist sdvDist meanSim sdvSim corr)]);
}
# Loop through them.
for my $genome0 (@$genomes) {
    # Are we resuming?
    if ($resume) {
        # Skip until we find the specified resume genome.
        $stats->Add(resumeSkip => 1);
        $nGroups++;
        if ($genome0 eq $resume) {
            $resume = '';
            print STDERR "Resuming after $genome0.\n" if $debug;
        }
    } else {
        # Get this genome's representative object and the genomes in its group.
        my $repGenome = $repDB->rep_object($genome0);
        $nGroups++;
        $stats->Add(groupIn => 1);
        my $repList = $repGenome->rep_list();
        my $name = $repGenome->name();
        my $listCount = scalar(@$repList) + 1;
        if ($listCount == 1) {
            print STDERR "Skipping singleton group $genome0: $name.\n" if $debug;
            $stats->Add(singletons => 1);
        } elsif ($listCount < $size) {
            print STDERR "Skipping small group ($listCount) $genome0: $name.\n" if $debug;
            $stats->Add(smalls => 1);
        } else {
            print STDERR "Processing $genome0 ($nGroups of $totalGroups) with $listCount members: $name.\n" if $debug;
            # This will count the number of times each role appears as a universal.
            my %counts;
            my @genomes = ($genome0, map { $_->[0] } @$repList);
            # Get the GEOs for the representatives.
            my $gHash = GEO->CreateFromPatric(\@genomes, %gOpts);
            # Loop through the members.
            for my $genome (sort keys %$gHash) {
                my $geo = $gHash->{$genome};
                my $roleH = $geo->roleCounts;
                for my $role (keys %$roleH) {
                    if ($roleH->{$role} == 1) {
                        $counts{$role}++;
                    }
                }
            }
            # Output the univeersal roles.
            my $found = 0;
            my $limit = $pct * $listCount / 100;
            print "$genome0\t$name\n";
            for my $role (sort keys %counts) {
                if ($counts{$role} >= $limit) {
                    print "\t$role\n";
                    $found++;
                    $stats->Add(foundRole => 1);
                }
            }
            print STDERR "$found universal roles in $genome0 group.\n" if $debug;
        }
    }
}
print STDERR "All done.\n" . $stats->Show() if $debug;