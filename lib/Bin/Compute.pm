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
    use IO::File;

=head1 Compute the Bins For a Community's Contigs.

This object contains the main method for converting contigs into bins. It takes as input a L<Bin::Score>
object containing scoring parameters, plus a list of bins, and an optional file containing scoring
vectors.

The fields of the object are as follows.

=over 4

=item logh

Open file handle for producing log and trace output.

=item stats

Statistics object for keeping statistics about a run.

=item score

L<Bin::Score> object for computing comparison scores between bins.

=back

=head2 Special Methods

    my $computer = Bin::Compute->new($score, %options);

Create a new bin computation object.

=over 4

=item score

A L<Bin::Score> object for comparing two bins.

=item options

A hash of optional parameters.

=over 8

=item logFile

An open file handle for the log/trace file, or the name of the log/trace file. The default is to write to STDOUT.

=back

=back

=cut

sub new {
    my ($class, $score, %options) = @_;
    # Handle the log file.
    my $logh;
    if (! $options{logFile}) {
        $logh = \*STDOUT;
    } elsif (! ref $options{logFile}) {
        $logh = IO::File->new(">$options{logFile}") || die "Could not open log: $!";
        $logh->autoflush(1);
    } else {
        $logh = $options{logFile};
    }
    # Create the object.
    my $retVal = {
        score => $score,
        logh => $logh,
        stats => Stats->new()
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Query Methods

=head3 stats

    my $stats = $helper->stats;

Return the L<Stats> object for tracking statistics.

=cut

sub stats {
    return $_[0]->{stats};
}


=head2 Public Methods

=head3 ProcessScores

    my $bins = $computer->ProcessScores($binList, $vectorFile);

Score comparisons between the specified bins to merge them into larger bins. We will first perform a comparison between
all the bins with found reference genomes. The non-zero comparisons will be sorted and used to combine the bins.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item vectorFile (optional)

If specified, a tab-delimited file containing comparison vectors for all the interesting contigs. For each pair of
interesting contigs, there will be a record containing (0) the first contig ID, (1) the second contig ID, (2) the
coverage score, (3) the tetranucleotide score, (4) the closest-reference-genome score, (4) the number of universal
roles not in common, and (5) the number of universal roles in common. If the file is not provided, the scores will
be computed from scratch.

=item RETURN

Returns a list of L<Bin> objects, sorted from most interesting to least, representing merged copies of the original
bins.

=back

=cut

sub ProcessScores {
    my ($self, $binList, $vectorFile) = @_;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Get the scoring object and the log file handle.
    my $score = $self->{score};
    my $logh = $self->{logh};
    # Get the number of bins.
    my $binCount = scalar @$binList;
    print $logh "$binCount input contigs specified for score processing.\n";
    $stats->Add(contigsRead => $binCount);
    # Select only the bins with reference genomes.
    my @selectedBins = grep { $_->hasHits } @$binList;
    my $filteredCount = scalar @selectedBins;
    my $totalScores = $filteredCount * ($filteredCount + 1) / 2;
    print $logh "$filteredCount useful contigs found in input. $totalScores scores required.\n";
    # This list contains 3-tuples for each pair of contigs containing the contig IDs and the comparison vector.
    my @scores;
    # This will count our progress.
    my $scoresComputed = 0;
    # Do we have a vector file?
    if ($vectorFile) {
        # Yes. Read from the file.
        print $logh "Reading score vectors from $vectorFile.\n";
        open(my $vh, "<", $vectorFile) || die "Could not open score vector file: $!";
        while (! eof $vh) {
            my $line = <$vh>;
            chomp $vh;
            my ($contigI, $contigJ, @vector) = split /\t/, $line;
            my $scoreValue = $score->ScoreV(\@vector);
            $stats->Add(pairsScored => 1);
            $scoresComputed++;
            if ($scoresComputed % 1000000 == 0) {
                print $logh "$scoresComputed of $totalScores scores computed.\n";
            }
            if ($scoreValue > 0) {
                $stats->Add(pairsKept => 1);
                push @scores, [$contigI, $contigJ, $scoreValue];
            }
        }
    } else {
        # No vector file. Compute the scores manually.
        my ($i, $j);
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
                    print $logh "$scoresComputed of $totalScores scores computed.\n";
                }
                $stats->Add(pairsScored => 1);
                if ($scoreValue > 0) {
                    $stats->Add(pairsKept => 1);
                    push @scores, [$contigI, $contigJ, $scoreValue];
                }
            }
        }
    }
    # We must loop through the contigs, comparing them. This hash tells us which bin
    # contains each contig.
    my %contig2Bin = map { $_->contig1 => $_ } @$binList;
    print $logh "Sorting " . scalar(@scores) . " scores.\n";
    @scores = sort { $b->[2] <=> $a->[2] } @scores;
    print $logh "Merging bins.\n";
    for my $scoreTuple (@scores) {
        # Get the bins.
        my ($contigI, $contigJ, $scoreValue) = @$scoreTuple;
        my $binI = $contig2Bin{$contigI};
        my $binJ = $contig2Bin{$contigJ};
        # Compute the score for these two bins.
        if (! $binI) {
            print $logh "WARNING: Missing bin for $contigI.\n";
            $scoreValue = 0;
            $stats->Add(missingBin => 1);
        } elsif (! $binJ) {
            print $logh "WARNING: Missing bin for $contigJ.\n";
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
    print $logh "Preparing for output.\n";
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
    print $logh scalar(@bins) . " bins output.\n";
    my @sortedBins = sort { Bin::cmp($a, $b) } @bins;
    # Return the statistics and the list of bins.
    return \@sortedBins;
}


1;