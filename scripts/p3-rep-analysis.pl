=head1 Analyze Taxonomic Characteristics of Representative-Genome Directories

    p3-rep-analysis.pl [options] repDir1 repDir2 repDir3 ... repDirN

This script will analyze the taxonomic composition of one or more representative-genome directories in an attempt to determine
how well they mirror known taxonomic classifications.

A specific set of genomes will be taken from the input.  For each genome, we will interrogate PATRIC for its major
taxonomic groupings.  Then we will run through the various representative-genome directories and perform the appropriate
taxonomic analysis.

Within a single representative-genome set, for each grouping, we want to know how many representative-genome groups contain genomes
from that rank.  For each representative-genome group we also want to know how many taxonomic groupings of each rank there are.

The supported ranks are the ones stored in the genome record-- B<kingdom>, B<phylum>, B<class>, C<order>, C<family>, C<genus>, C<species>.

=head2 Parameters

The positional parameters are the names of the representative genome directories.  Each will be processed individually, and a
tabular summary produced on the standard output.

The standard input should contain the incoming genome IDs and can be overridden using the options in L<P3Utils/ih_options>.

Additional command-line options are those given in L<P3Utils/col_options> (to specify the genome ID column in the input file) plus
the following options.

=over 4

=item verbose

Progress messages will be produced on STDERR.

=item save

If specified, then the group assignments computed by this program will be saved in a C<rep_db.new.tbl> file in each representative
server directory.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use RoleParse;
use Math::Round;

use constant RANKS => [qw(kingdom phylum class order family genus species)];

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir1 repDir2 repDir3 ... repDirN', P3Utils::ih_options(),
        P3Utils::col_options(),
        ['verbose|debug|v', 'show progress on STDERR'],
        ['save', 'save group assignments']
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the options.
my $debug = $opt->verbose;
my $save = $opt->save;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
print STDERR "Reading genome input.\n" if $debug;
# Compute the return columns for the genomes.
my @cols = (qw(genome_id genome_name), @{RANKS()});
# Compute the parameters for the seed protein retrieval.
my @filter = (['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']);
my @cols2 = qw(genome_id genome_name product aa_sequence);
# Save the checksum for the seed role.
my $roleCheck = "WCzieTC/aZ6262l19bwqgw";
# This is a two-level hash, keyed by genome ID and then taxonomic rank.  The value is the grouping name for the given genome
# at the given rank.
my %genomes;
# This is a single-level hash, keyed by genome ID.  For each genome it returns the seed protein value.
my %prots;
# Loop through the input.
my $batchCount = 0;
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    $batchCount++;
    print STDERR "Processing genome batch $batchCount.\n" if $debug;
    # First we get the taxonomic information for this batch.
    my $genomeList = [map { $_->[0] } @$couplets];
    my $resultList = P3Utils::get_data_keyed($p3, genome => [], \@cols, $genomeList);
    for my $result (@$resultList) {
        my ($genomeID, $name, @taxons) = @$result;
        for (my $i = 0; $i < scalar @taxons; $i++) {
            my $rank = RANKS->[$i];
            my $taxon = $taxons[$i];
            if (! $taxon) {
                print STDERR "Missing $rank for $genomeID: $name.\n" if $debug;
                $taxon = '<unknown>';
            }
            $genomes{$genomeID}{$rank} = $taxon;
        }
    }
    # Now we get the seed proteins for this batch.
    $resultList = P3Utils::get_data_keyed($p3, 'feature', \@filter, \@cols2, $genomeList, 'genome_id');
    # The resultList entries are in the form [$genome, $name, $function, $prot]. Get the longest
    # protein for each genome. We also warn about obviously bad genomes.
    for my $result (@$resultList) {
        my ($genome, $name, $function, $prot) = @$result;
        # Check the protein.
        my $check = RoleParse::Checksum($function // '');
        if (! $prot) {
            print STDERR "WARNING: $genome $name has no identifying protein.\n" if $debug;
        } elsif ($check eq $roleCheck) {
            # Here the function matched and is correct.
            $prots{$genome} = $prot;
        }
    }
}
# Now we know the taxonomy of each genome at each rank level.  For each representative-genome directory, we want to know the
# average number of taxonomic groups per representative group, and the average number of representative groups per taxonomic group.
# We don't use the precomputed representations; rather, we find the absolutely closest representative for each genome using the
# seed protein.  Note also that when we compute the means, we ignore the singleton groups.
#
# This is a two-level hash, keyed by repdir name and then rank.  The value is the number of rep-group/tax-group pairs in that
# directory for that rank.  There is an extra rank for genomes themselves, "member".
my %pairings;
# This is a two-level hash, also keyed by repdir name and then rank.  The value is the number of tax groups in that directory
# for that rank.  There is an extra rank for genomes themselves, "member".
my %taxGroups;
# This is a single-level hash, keyed by repdir.  The value is the number of rep groups in that directory.
my %repGroups;
# This is a single-level hash, keyed by repdir.  The value is the number of singleton rep groups in that directory.
my %singletons;
# This is a single-level hash, keyed by repdir.  The value is the member count of the largest rep group in that directory.
my %largest;
# Now we can loop through the directories.
for my $repDir (@ARGV) {
    print STDERR "Processing representative directory $repDir.\n" if $debug;
    my $repDB = RepGenomeDb->new_from_dir($repDir, unconnected => 1);
    my $score = $repDB->score();
    print STDERR "Similarity score for this directory is $score.\n" if $debug;
    # If we are saving, open the save file.
    my $oh;
    if ($save) {
        open($oh, ">$repDir/rep_db.new.tbl") || die "Could not open output rep_db: $!";
    }
    # This hash is keyed by rank and then pairing (name\tgroup).  It enables us to count the pairings.
    my %rPairings;
    # This hash is keyed by group and counts the genomes in each group.
    my %groups;
    # This hash is keyed by rank and then taxon name.  It enables us to count the taxa per rank.
    my %taxa;
    # Loop through all the genomes, comparing them to this representative genome set.
    my $gCount = 0;
    for my $genome (keys %prots) {
        my ($repID, $score2) = $repDB->find_rep($prots{$genome});
        $gCount++;
        if ($score2 >= $score) {
            # Here the genome is in a group.
            $groups{$repID}++;
            # Update the pair lists.
            for my $rank (@{RANKS()}) {
                my $taxon = $genomes{$genome}{$rank};
                $rPairings{$rank}{"$taxon\t$repID"} = 1;
                $taxa{$rank}{$taxon} = 1;
            }
            # Save if necessary.
            if ($oh && $genome ne $repID) {
                P3Utils::print_cols([$genome, $repID, $score2], oh => $oh);
            }
        }
        print STDERR "$gCount genomes processed.\n" if $debug && ($gCount % 100 == 0);
    }
    print STDERR "Computing totals.\n" if $debug;
    # Now find the largest group and the number of singletons.
    my ($largest, $largestCount, $singletons, $groupCount, $memCount) = ('', 0, 0, 0);
    for my $repID (keys %groups) {
        my $size = $groups{$repID};
        $memCount += $size;
        if ($size == 1) {
            $singletons++;
        }
        if ($size > $largestCount) {
            $largestCount = $size;
            $largest = $repID;
        }
        $groupCount++;
    }
    print STDERR "Largest group is $largest ($largestCount members).  $singletons singleton groups out of $groupCount.\n" if $debug;
    $largest{$repDir} = $largestCount;
    $singletons{$repDir} = $singletons;
    $repGroups{$repDir} = $groupCount;
    my (%gPairings, %gTaxa);
    $gPairings{member} = $memCount;
    $gTaxa{member} = $memCount;
    for my $rank (@{RANKS()}) {
        $gPairings{$rank} = scalar keys %{$rPairings{$rank}};
        $gTaxa{$rank} = scalar keys %{$taxa{$rank}};
    }
    $pairings{$repDir} = \%gPairings;
    $taxGroups{$repDir} = \%gTaxa;
}
print STDERR "Producing reports.\n" if $debug;
# Now we produce the reports.  First, we need a header line based on the the repDir names.
P3Utils::print_cols(['Statistic', @ARGV]);
# First row is number of groups.
P3Utils::print_cols(['groups', map { $repGroups{$_} } @ARGV]);
# Second row is number of singletons.
P3Utils::print_cols(['singletons', map { $singletons{$_} } @ARGV]);
# Third row is largest group size.
P3Utils::print_cols(['largest', map { $largest{$_} } @ARGV]);
# Now the mean number of groups per rank.
for my $rank (@{RANKS()}) {
    my %means;
    for my $repDir (@ARGV) {
        my $denom = $repGroups{$repDir} - $singletons{$repDir};
        my $mean = 1;
        if ($denom > 0) {
            $mean = Math::Round::nearest(0.1, ($pairings{$repDir}{$rank} - $singletons{$repDir}) / $denom);
        }
        $means{$repDir} = $mean;
    }
    P3Utils::print_cols(["mean $rank/group", map { $means{$_} } @ARGV]);
}
# Now the mean number of ranks per group.
for my $rank (@{RANKS()}) {
    my %means;
    for my $repDir (@ARGV) {
        my $denom = $taxGroups{$repDir}{$rank};
        my $mean = 1;
        if ($denom > 0) {
            $mean = $pairings{$repDir}{$rank} / $denom;
        }
        $means{$repDir} = Math::Round::nearest(0.1, $mean);
    }
    P3Utils::print_cols(["mean groups/$rank", map { $means{$_} } @ARGV]);
}
print STDERR "All done.\n" if $debug;