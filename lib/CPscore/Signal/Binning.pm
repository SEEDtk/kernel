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


package CPscore::Signal::Binning;

    use strict;
    use warnings;
    use base qw(CPscore::Signal);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on the
precise ordering of the signals in similarity signal vectors.

=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    my $retVal = $self->SUPER::type() . 'bin';
    return $retVal;
}


=head3 vector_compare

    my $score = $scoring->vector_compare($contig1, $cv1, $contig2, $cv2);

Compute the similarity score for a pair of contigs from their similarity vectors. Each similarity
vector contains the scores computed by L<CPscore::Vector/update_vector_hash> in order by the relevant reference
genome ID and formatted by L<CPscore/store_vector>. The score is C<1> if the two vectors have the same
ordering and C<0> otherwise.

=over 4

=item contig1

ID of the first contig.

=item cv1

Similarity vector for the first contig.

=item contig2

ID of the second contig.

=item cv2

Similarity vector for the second contig.

=item RETURN

Returns the similarity score for the two contigs.

=back

=cut

sub vector_compare {
    # Get the parameters.
    my ($self, $contig1, $cv1, $contig2, $cv2) = @_;
    # We need to find the coordinates in common. We put them into the following vector in the form
    # [score1, score2].
    my @scores;
    my $n = scalar @$cv1;
    for (my $i = 0; $i < $n; $i++) {
        push @scores, [$cv1->[$i], $cv2->[$i]];
    }
    my @sorted = sort { $b->[0] <=> $a->[0] } @scores;
    # The vector is now sorted by the CV1 scores, in order from largest to smallest.
    # We will now compare every two positions. If a subsequent position has a greater
    # CV2 score than a particular position, we have a failure and return 0.
    my $retVal = 1;
    $n = scalar(@sorted) - 1;
    # Loop through the pairs, searchign for a mismatch.
    for (my $i = 0; $i < $n && $retVal; $i++) {
        my $sortedI = $sorted[$i][1];
        my $sortedJ = $sorted[$i+1][1];
        if ($sortedJ > $sortedI) {
            $retVal = 0;
        }
    }
    # Return the tally.
    return $retVal;
}


1;