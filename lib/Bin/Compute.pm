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


package Bin::Compute;

    use strict;
    use warnings;
    use Bin;
    use Stats;

=head1 Compute the Bins For a Community's Contigs.

This package contains the main method for converting contigs into bins. It takes as input a L<Bin::Score>
object containing scoring parameters, plus an open input file handle for the contig data.

=head2 Public Methods

=head3 ScoreVectors

    my ($stats, $vectors) = Bin::Compute::ScoreVectors($binList)'

Compute the scoring vectors for a set of contigs. For each likely pair of contigs, we compute the vector
of scoring components. These vectors can be used to compute actual scores.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item RETURN

Returns a two-element list consisting of (0) a L<Stats> object containing statistics about the scoring and
(2) a list of 3-tuples, each containing [0] a scoring vector, [1] the first contig ID, and [2] the second
contig ID.

=back

=cut

sub ScoreVectors {
    my ($binList) = @_;
    # Create the statistics object.
    my $stats = Stats->new();
    my $binCount = scalar @$binList;
    $stats->Add(contigsRead => $binCount);
    # Select only the bins with reference genomes.
    my @selectedBins = grep { $_->hasHits } @$binList;
    my $filteredCount = scalar @selectedBins;
    my $totalScores = $filteredCount * ($filteredCount + 1) / 2;
    print "$filteredCount useful contigs found in input. $totalScores scores required.\n";
    # This list contains 3-tuples for each pair of contigs containing the contig IDs and the comparison vector.
    my @scores;
    my ($i, $j);
    my $scoresComputed = 0;
    # Now the nested loop.
    for ($i = 0; $i < $filteredCount; $i++) {
        my $binI = $selectedBins[$i];
        my $contigI = $binI->contig1;
        for ($j = $i + 1; $j < $filteredCount; $j++) {
            my $binJ = $selectedBins[$j];
            my $contigJ = $binJ->contig1;
            my $scoreVector = Bin::Score::Vector($binI, $binJ);
            $scoresComputed++;
            if ($scoresComputed % 1000000 == 0) {
                print "$scoresComputed of $totalScores scores computed.\n";
            }
            $stats->Add(pairsScored => 1);
            push @scores, [$scoreVector, $contigI, $contigJ];
        }
    }
    # Return the statistics and the score vector list.
    return ($stats, \@scores);
}


=head3 ProcessScores

    my ($stats, $bins) = Bin::Compute::ProcessScores($binList, $score, $scoreVectors);

Process a set of contigs and form them into bins.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item score

A L<Bin::Score> object containing the scoring parameters.

=item scoreVectors (optional)

Reference to a list of 3-tuples produced by L<ScoreVectors>, each 3-tuple consisting of (0) a scoring vector,
(1) the first contig ID, and (2) the second contig ID.

=item RETURN

Returns a two-element list consisting of (0) a L<Stats> object containing statistics for the operation and (1) a
reference to a list of the output bins. If it is not supplied it will be computed internally.

=cut

sub ProcessScores {
    my ($binList, $score, $scoreVectors) = @_;
    # Create the statistics object.
    my $stats = Stats->new();
    # This list will contains 3-tuples for each comparable pair of contigs containing the contig IDs and the comparison score.
    my @scores;
    # We have two possible scenarios: scoring vectors were input, or we compute scores on the fly. The scoring vectors are very
    # memory-intensive.
    if ($scoreVectors) {
        my $incomingScores = scalar @$scoreVectors;
        print "Converting $incomingScores score components to final scores.\n";
        for my $scoreVector (@$scoreVectors) {
            my ($vector, $contigI, $contigJ) = @$scoreVector;
            my $scoreValue = $score->ScoreV($vector);
            $stats->Add(scoresComputed => 1);
            if ($scoreValue > 0) {
                $stats->Add(pairsKept => 1);
                push @scores, [$contigI, $contigJ, $scoreValue];
            }
        }
    } else {
        # Select only the bins with reference genomes.
        my @selectedBins = grep { $_->hasHits } @$binList;
        my $filteredCount = scalar @selectedBins;
        my $totalScores = $filteredCount * ($filteredCount + 1) / 2;
        print "$filteredCount useful contigs found in input. $totalScores scores required.\n";
        my ($i, $j);
        my $scoresComputed = 0;
        # Now the nested loop.
        for ($i = 0; $i < $filteredCount; $i++) {
            my $binI = $selectedBins[$i];
            my $contigI = $binI->contig1;
            for ($j = $i + 1; $j < $filteredCount; $j++) {
                my $binJ = $selectedBins[$j];
                my $contigJ = $binJ->contig1;
                my $scoreValue = $score->Score($binI, $binJ);
                $scoresComputed++;
                if ($scoresComputed % 1000000 == 0) {
                    print "$scoresComputed of $totalScores scores computed.\n";
                }
                $stats->Add(pairsScored => 1);
                if ($scoreValue > 0) {
                    push @scores, [$contigI, $contigJ, $scoreValue];
                }
            }
        }
    }
    # We must loop through the contigs, comparing them. This hash tells us which bin
    # contains each contig.
    my %contig2Bin = map { $_->contig1 => $_ } @$binList;
    print "Sorting " . scalar(@scores) . " scores.\n";
    @scores = sort { $b->[2] <=> $a->[2] } @scores;
    print "Merging bins.\n";
    for my $scoreTuple (@scores) {
        # Get the bins.
        my ($contigI, $contigJ, $scoreValue) = @$scoreTuple;
        my $binI = $contig2Bin{$contigI};
        my $binJ = $contig2Bin{$contigJ};
        # Compute the score for these two bins.
        if (! $binI) {
            print "WARNING: Missing bin for $contigI.\n";
            $scoreValue = 0;
            $stats->Add(missingBin => 1);
        } elsif (! $binJ) {
            print "WARNING: Missing bin for $contigJ.\n";
            $scoreValue = 0;
            $stats->Add(missingBin => 1);
        } elsif ($binI->contig1 eq $binJ->contig1) {
            # Here the contigs are already in the same bin.
            $scoreValue = 0;
            $stats->Add(binsAlreadyMerged => 1);
        } elsif ($binI->contigCount > 1 || $binJ->contigCount > 1) {
            # The bins are not singletons. We have to re-compute the score.
            $stats->Add(complexBinMatch => 1);
            $scoreValue = $score->Score($binI, $binJ);
        }
        # If the score is nonzero, we can merge them.
        if ($scoreValue > 0) {
            $stats->Add(binsMerged => 1);
            # The bins will be merged into bin I. First, we must update the contig-to-bin map.
            for my $contigJX ($binJ->contigs) {
                $contig2Bin{$contigJX} = $binI;
            }
            # Now merge the bins.
            $binI->Merge($binJ);
        }
    }
    # Get the final list of bins.
    print "Preparing for output.\n";
    my %found;
    my @bins;
    for my $contigX (keys %contig2Bin) {
        my $bin = $contig2Bin{$contigX};
        my $binID = $bin->contig1;
        if (! $found{$binID}) {
            push @bins, $bin;
            $found{$binID} = 1;
            $stats->Add(outputBin => 1);
        }
    }
    print scalar(@bins) . " bins output.\n";
    my @sortedBins = sort { Bin::cmp($a, $b) } @bins;
    # Return the statistics and the list of bins.
    return ($stats, \@sortedBins);
}

1;