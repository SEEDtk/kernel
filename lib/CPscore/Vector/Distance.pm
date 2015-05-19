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


package CPscore::Vector::Distance;

    use strict;
    use warnings;
    use base qw(CPscore::Vector);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on the
relative distance of the signals in similarity signal vectors. The distance in this case is not a true metric
distance. Instead, the vector is scored positively for matching nonzero coordinates and negatively for mismatches.
If both coordinates are zero, they are ignored. Otherwise, the absolute value of the difference is taken
and the score is incremented by 1/(1+I<x>), where I<x> is the percent difference.

=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    my $retVal = $self->SUPER::type() . 'dist';
    return $retVal;
}


=head3 vector_compare

    my $score = $scoring->vector_compare($contig1, $cv1, $contig2, $cv2);

Compute the similarity score for a pair of contigs from their similarity vectors. Each similarity
vector contains the scores computed by L<CPscore::Vector/update_vector_hash> in order by the relevant reference
genome ID and formatted by L<CPscore/store_vector>.

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
    # Loop through the two vectors.
    my $retVal = 0;
    my $n = scalar @$cv1;
    for (my $i = 0; $i < $n; $i++) {
        my $v1 = $cv1->[$i];
        my $v2 = $cv2->[$i];
        if ($v1 || $v2) {
            my $x = (($v1 > $v2) ? ($v1 - $v2)/$v1 : ($v2 - $v1)/$v2);
            $retVal += 1 / (1 + $x);
        }
    }
    # Return the tally.
    return $retVal;
}


1;