=head1 Compute the Universal Role Distance Between Represented Genomes

    p3-rep-role-distance.pl [options] repDir

This program will analyze the represented-genome groups in a representative-genome directory.  For each represented set, it
will return the mean, median, minimum, maximum, and standard deviation of the universal-role distance (which varies from 0 to 100),
and the mean and standard deviation of the similarity scores.  It will also return the number of genomes in the group.  The distances
are all computed from the representative genome.

=head2 Parameters

The positional parameter is the name of the representative-genome directory.

Additional command-line options are those in L<EvalCon/role_options> plus the following.

=over 4

=item verbose

Display progress messages on STDERR.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use Statistics::Basic qw(:all);
use EvalCon;
use Stats;
use GEO;


# Get the command-line options.
my $opt = P3Utils::script_opts('repDir', EvalCon::role_options(),
        ['verbose|debug|v', 'display progress messages on STDERR'],
        );
# Get the options.
my $debug = $opt->verbose;
my $logH = ($debug ? \*STDERR : undef);
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
# Set up the formatter.
my $formatter = Number::Format->new();
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
P3Utils::print_cols([qw(group name size minDist meanDist maxDist sdvDist meanSim sdvSim corr)]);
# Loop through them.
for my $genome0 (@$genomes) {
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
    } else {
        print STDERR "Processing $genome0 ($nGroups of $totalGroups) with $listCount members: $name.\n" if $debug;
        my %sims = map { $_->[0] => $_->[1] } @$repList;
        # Get the GEOs for the representatives.
        my $gHash = GEO->CreateFromPatric([$genome0, keys %sims], %gOpts);
        # Get the basic GEO.
        my $geo0 = $gHash->{$genome0};
        if (! $geo0) {
            print STDERR "Genome error on $genome0.\n" if $debug;
            $stats->Add(groupSkipped => 1);
        } else {
            # Loop through the members.
            my (@distances, @sims);
            my ($minDist, $maxDist) = (100, 0);
            for my $genome (sort keys %$gHash) {
                if ($genome ne $genome0) {
                    my $geo = $gHash->{$genome};
                    my $distance = $geo0->uni_similarity($geo);
                    push @distances, $distance;
                    push @sims, $sims{$genome};
                    $stats->Add(distanceComputed => 1);
                    $minDist = $distance if $distance < $minDist;
                    $maxDist = $distance if $distance > $maxDist;
                }
            }
            print STDERR "Calculating statistics for $genome0.\n" if $debug;
            # Now compute the stats and write them.
            my $dV = vector(@distances);
            my ($mean, $std) = (mean($dV), stddev($dV));
            my $sV = vector(@sims);
            my ($meanSim, $stdSim) = (mean($sV), stddev($sV));
            my $corr = correlation($dV, $sV);
            P3Utils::print_cols([$genome0, $name, $listCount,
                    map { $formatter->format_number($_, 2, 1) } $minDist, $mean, $maxDist, $std, $meanSim, $stdSim, $corr]);
        }
    }
}
print STDERR "All done.\n" . $stats->Show() if $debug;