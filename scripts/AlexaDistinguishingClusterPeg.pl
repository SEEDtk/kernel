=head1 Alexa Job to Find a Central Signature Cluster Peg

    AlexaDistinguishingClusterPeg [options]

This script takes as input two genome sets. Its purpose is to find a candidate for compare regions analysis that will
reveal information about a signature cluster that distinguishes between two genome sets. We first run multiple
signature-cluster analyses on subsets of the incoming genome sets. Random genomes will be chosen from the two
sets and signature families computed for each set pair, followed by clusters of the families. After running the
set of iterations, we compute the number of times each pair of signature families occurs in a signature cluster.

It is important to understand that these pairs represent families that occur in the first input genome set but
not the second.

This script produces an output table and an output set. The table has three columns-- count, family1.family_id, and
family2.family_id-- representing the number of times a family pair occurred and the two famiies in the pair. This
table is sorted in descending order by count.

The set contains the pegs of interest. The default is one peg, but multiple pegs can be requested.  

=head2 Parameters

The positional parameters are the two input genome sets, and the output table and set.  Only the set/table name
is specified, not the file name (e.g. C<F> instead of C<F.tbl>). A single output name is specified, but both an
output set and an output table will be produced.

The typical command-line parameters required by Alexa jobs are supported. These are documented in the L<Job> object.

The following command-line parameters are supported by this script.

=over 4

=item size

Sample size to use when picking items from the sets.

=item iterations

Number of iterations to run.

=item pegs

Number of output features to select.

=item minIn

Minimum fraction of genomes in set 1 that must contain a signature family.

=item maxOut

Maximum fraction of genomes in set 2 that may contain a signature family.

=back

=cut

use strict;
use Job_Config;
use Job;
use P3Signatures;
use SeedAware;
use P3DataAPI;
use P3Utils;

# Parse the command line and create the job object.
my $jobObject = Job->new('set1 set2 result', 
        ['size=i', 'number of genomes to use per set in each iteration', { default => 20 }],
        ['iterations=i', 'number of iterations to run', { default => 10 }],
        ['pegs=i', 'number of pegs to return', { default => 1 }],
        ['minIn=f', 'minimum fraction of set 1 genomes that must contain a signature family', { default => 1 }],
        ['maxOut=f', 'maximum fraction of set 2 genomes that may contain a signature family', { default => 0 }]
        );
# These variables are needed outside the EVAL block.
my $count = 0;
my $iter = 0;
# Protect against errors.
eval {
    # Get the command-line options.
    my $size = $jobObject->opt->size;
    my $iterations = $jobObject->opt->iterations;
    my $pegs = $jobObject->opt->pegs;
    my $min_in = $jobObject->opt->minin;
    my $max_out = $jobObject->opt->maxout;
    # Connect to PATRIC.
    my $p3 = P3DataAPI->new();
    # Get the working directory.
    my $workDir = $jobObject->workDir;
    # Get sets and table.
    my ($set1, $set2, $tableO) = @ARGV;
    my $setO = $tableO;
    die "Insufficient parameters for job." if ! $tableO;
    # Read the input sets.
    my $genomes1L = $jobObject->ReadSet($set1);
    my $genomes2L = $jobObject->ReadSet($set2);
    # Compute the shortest set.
    my $setLength = scalar @$genomes1L;
    if (scalar(@$genomes2L) < $setLength) {
        $setLength = scalar(@$genomes2L)
    }
    # Insure our size isn't too big.
    $size = $setLength if ($size > $setLength);
    # Position at the end of the set to force an initial shuffle.
    my $position = $setLength;
    # This hash will track the number of times each family pair appears in the output.
    my %pairCounts;
    # This hash contains the clusters inside each family.
    my %familyClusters;
    # Loop through the iterations.
    for ($iter = 1; $iter <= $iterations; $iter++) {
        # Do we need to reshuffle the input sets?
        my $end = $position + $size;
        if ($end > $setLength) {
            # Yes. Shuffle them now. Note we choose random items from the whole set, but only shuffle
            # into the length we're using.
            for my $set ($genomes1L, $genomes2L) {
                my $total = scalar @$set;
                for (my $i = 0; $i < $setLength; $i++) {
                    my $j = int(rand($total));
                    ($set->[$i], $set->[$j]) = ($set->[$j], $set->[$i])
                }
            }
            # Reposition at the beginning.
            $position = 0;
            $end = $size;
        }
        # Get the current slice of each set.
        my @subset1 = @{$genomes1L}[$position .. ($end - 1)];
        my @subset2 = @{$genomes2L}[$position .. ($end - 1)];
        # Compute the signature families.
        my $familyHash = P3Signatures::Process(\@subset1, \@subset2, $min_in, $max_out, $jobObject, "Iteration $iter.");
        # Compute the fields to select.
        my @selectList = qw(patric_id accession start end strand product);
        # Now we create the peg-info file.
        $jobObject->Progress("Creating peg information file for iteration $iter.");
        my $pegInfoFile = "$workDir/peginfo.tbl";
        open(my $oh, ">$pegInfoFile") || die "Could not open peg info file: $!";
        print $oh join("\t", 'family.family_id', map { "feature.$_" } @selectList) . "\n";
        # Compute the genome filter.
        my $gFilter = ['in', 'genome_id', '(' . join(',', @subset1) . ')'];
        # Loop through the hash, creating couplets and submitting them in batches for feature data.
        my @couplets;
        for my $fam (keys %$familyHash) {
            if (scalar @couplets > 100) {
                my $batch = P3Utils::get_data_batch($p3, feature => [$gFilter], \@selectList, \@couplets, 'plfam_id');
                for my $line (@$batch) {
                    P3Utils::print_cols($line, $oh);
                }
                @couplets = ();
            } else {
                push @couplets, [$fam, [$fam]];
            }
        }
        my $batch = P3Utils::get_data_batch($p3, feature => [$gFilter], \@selectList, \@couplets, 'plfam_id');
        for my $line (@$batch) {
            P3Utils::print_cols($line, $oh);
        }
        close $oh;
        # Release the family hash memory.
        undef $familyHash;
        # Compute the signature clusters.
        $jobObject->Progress("Computing signature clusters for iteration $iter.");
        my $clusterFile = "$workDir/clusters.tbl";
        my @clusterLines = `p3-signature-clusters --input=$pegInfoFile`;
        # The clusters list contains entries of the form "familyID\tpeg\tfunction" for each peg in a cluster, with the
        # clusters themselves separated by double slashes. We will accumulate the current cluster in this list, in the
        # form [family-id,peg] for each peg.
        $jobObject->Progress("Analyzing signature clusters for iteration $iter.");
        my @pegCluster;
        my @familyCluster;
        for my $line (@clusterLines) {
            chomp $line;
            if ($line ne '//') {
                my ($family, $peg) = split /\t/, $line;
                push @pegCluster, $peg;
                push @familyCluster, $family;
            } else {
                # Here we have a whole cluster. Count the pairs.
                my $n = (scalar @pegCluster);
                # Get a safe copy of the peg list. We will be storing this in the family cluster hash.
                my $pegList = [@pegCluster];
                # Loop through the families, processing pairs.
                for (my $i = 0; $i < $n; $i++) {
                    my $family = $familyCluster[$i];
                    # Associate the list of pegs with the  family.
                    push @{$familyClusters{$family}}, $pegList;
                    # Process the pairs.
                    for (my $j = $i + 1; $j < $n; $j++) {
                        my $family2 = $familyCluster[$j];
                        if ($family2 ne $family) {
                            my $pair = join("\t", sort($family, $family2));
                            $pairCounts{$pair}++;
                        }
                    }
                }
                # Restart the cluster;
                @familyCluster = ();
                @pegCluster = ();
            }
        }
        # Push the position forward.
        $position = $end;
        # Count this iteration.
        $count++;
    }
    # Sort the pairs.
    $jobObject->Progress("Sorting family pairs.");
    my @pairs = sort { $pairCounts{$b} <=> $pairCounts{$a} } keys %pairCounts;
    # Output the table.
    $jobObject->StoreTable($tableO, "cluster pairs from signature families distinguishing $set1 from $set2",
        ['family1', 'family2', 'count'], \%pairCounts, \@pairs);
    # Loop through the pairs from which we want to output pegs. We use a splice to chop off the pairs that are
    # not of interest.
    $jobObject->Progress("Searching for pegs of interest.");
    my @outPegs;
    splice @pairs, $pegs;
    for my $pair (@pairs) {
        # Get the first family of the pair.
        my ($family) = split /\t/, $pair;
        # Get all of its clusters.
        my $clusterList = $familyClusters{$family};
        # Extract the longest cluster.
        my $cluster = [];
        my $clusterL = 0;
        for my $clusterI (@$clusterList) {
            my $clusterIL = scalar @$clusterI;
            if ($clusterIL > $clusterL) {
                $cluster = $clusterI;
                $clusterL = $clusterIL;
            }
        }
        # Get the middle peg.
        my $peg = $cluster->[int($clusterL/2)];
        push @outPegs, $peg;
    }
    # Output the pegs found.
    $jobObject->StoreSet($setO, "interesting features from signature families distinguishing $set1 from $set2", \@outPegs);
};
if ($@) {
    $jobObject->Fail("ERROR during iteration $iter: $@");
} else {
    $jobObject->Finish("$count iterations completed.");
}