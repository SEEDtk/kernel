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


package CPscore::Compare::Distance;

    use strict;
    use warnings;
    use base qw(CPscore::Compare);

=head1 Community Pipeline Vector Scoring -- Distance Comparison

This is a community pipeline scoring object that computes the similarity between two sample contigs based on the
relative distance of the signals in similarity signal vectors. The distance in this case is not a true metric
distance. Instead, the vector is scored according to the degree of match between coordinates.
If both coordinates are zero, they are ignored. Otherwise, the absolute value of the difference is taken
and the score is incremented by 1/(1+I<x>), where I<x> is the percent difference.


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