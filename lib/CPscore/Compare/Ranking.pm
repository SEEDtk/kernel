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


package CPscore::Compare::Ranking;

    use strict;
    use warnings;
    use base qw(CPscore::Compare);

=head1 Community Pipeline Vector Scoring -- Ranking Comparison

This is a community pipeline scoring object that computes the similarity between two sample contigs based on the
relative strength of the signals in similarity signal vectors.



=head2 Virtual Methods

=head3 compare

    my $sim = $compare->compare($contig1, $contig2);

Return the similarity score for a pair of contigs.

=over 4

=item contig1

L<SampleContig> object for the first contig.

=item contig2

L<SampleContig> object for the second contig.

=item RETURN

Returns the similarity score for the two contigs.

=back

=cut

sub compare {
    # Get the parameters.
    my ($self, $contig1, $contig2) = @_;
    # Get the scoring vectors.
    my $cv1 = $contig1->vector;
    my $cv2 = $contig2->vector;
    # We need to find the coordinates in common. We put them into the following vector in the form
    # [score1, score2].
    my @scores;
    my $n = scalar @$cv1;
    for (my $i = 0; $i < $n; $i++) {
        if ($cv1->[$i] && $cv2->[$i]) {
            push @scores, [$cv1->[$i], $cv2->[$i]];
        }
    }
    my @sorted = sort { $b->[0] <=> $a->[0] } @scores;
    # The vector is now sorted by the CV1 scores, in order from largest to smallest.
    # We will now compare every two positions. If a subsequent position has a greater
    # CV2 score than a particular position, we count that as an ordering mismatch.
    # Otherwise it is an ordering match.
    my $retVal = 0;
    $n = scalar @sorted;
    if ($n > 0) {
        # If we have at least one matching coordinate, start the tally at 1.
        $retVal = 1;
        # Loop through the pairs, tallying matches.
        for (my $i = 0; $i < $n; $i++) {
            my $sortedI = $sorted[$i][1];
            for (my $j = $i + 1; $j < $n; $j++) {
                my $sortedJ = $sorted[$j][1];
                if ($sortedJ > $sortedI) {
                    $retVal -= 2;
                } else {
                    $retVal++;
                }
            }
        }
    }
    # Return the tally.
    return $retVal;
}


1;