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

=head3 Process

    my ($stats, $bins) = Bin::Compute::Process($binList, $score);

Process a set of contigs and form them into bins.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item score

A L<Bin::Score> object containing the scoring parameters.

=item RETURN

Returns a two-element list consisting of (0) a L<Stats> object containing statistics for the operation and (1) a
reference to a list of the output bins.

=cut

sub Process {
    my ($binList, $score) = @_;
    # Create the statistics object.
    my $stats = Stats->new();
    my $binCount = scalar @$binList;
    $stats->Add(contigsRead => $binCount);
    print "$binCount contigs found in input.\n";
    # We must loop through the contigs, comparing them. This hash tells us which bin
    # contains each contig.
    my %contig2Bin = map { $_->contig1 => $_ } @$binList;
    # This list contains 3-tuples for each pair of contigs containing the contig IDs and the comparison score.
    my @scores;
    my ($i, $j);
    # Now the nested loop.
    for ($i = 0; $i < $binCount; $i++) {
        my $binI = $binList->[$i];
        my $contigI = $binI->contig1;
        print "Processing scores for $contigI.\n";
        for ($j = $i + 1; $j < $binCount; $j++) {
            my $binJ = $binList->[$j];
            my $contigJ = $binJ->contig1;
            my $scoreValue = $score->($binI, $binJ);
            $stats->Add(pairsScored => 1);
            if ($scoreValue > 0) {
                $stats->Add(pairsKept => 1);
                push @scores, [$contigI, $contigJ, $scoreValue];
            }
        }
    }
    print "Sorting scores.\n";
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
    # Return the statistics and the list of bins.
    return ($stats, \@bins);
}

1;