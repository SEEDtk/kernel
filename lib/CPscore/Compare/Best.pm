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


package CPscore::Compare::Best;

    use strict;
    use warnings;
    use base qw(CPscore::Compare);

=head1 Community Pipeline Scoring -- Best Coordinate Comparison

This is a community pipeline scoring object that computes the similarity between two sample contigs based on a
conformity of the top-ranked signal strengths. The similarity score is I<1/(1+x)>, where I<x> is the percent
difference between the top coordinates in each vector. If the top coordinates are different, the
similarity score is 0.


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
    my $retVal = 0;
    # Check the two vectors for a coordinate match.
    my $i1 = $contig1->bestCoord;
    my $i2 = $contig2->bestCoord;
    if ($i1 == $i2) {
        # The coordinates match, so compare them.
        my $v1 = $contig1->vector($i1);
        my $v2 = $contig2->vector($i2);
        if ($v1 > $v2) {
            $retVal = 1 / (1 + ($v1 - $v2)/$v1);
        } else {
            $retVal = 1 / (1 + ($v2 - $v1)/$v2);
        }
    }
    # Return the score.
    return $retVal;
}


1;