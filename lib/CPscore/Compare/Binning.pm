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


package CPscore::Compare::Binning;

    use strict;
    use warnings;
    use base qw(CPscore::Compare);

=head1 Community Pipeline Scoring -- Bin Comparison

This is a community pipeline scoring object that computes the similarity between two sample contigs based on
whether they belong to the same ordering bin. The score is C<1> if the coordinates in the vectors are in the
same relative order-- that is, if the highest coordinate is at the same position in both, as is the second
highest, and so forth. Otherwise, it is C<0>.

In addition to the fields in the base class, this object contains the following.

=over 4

=item topSize

The maximum number of vector positions to consider when determining the bin. The default is 0, indicating all
positions are considered.

=back

=head2 Special Methods

=head3 new

    my $compare = CPscore::Compare->new(%options);

Create a new, blank comparison object.

=over 4

=item options

Hash containing options for the comparison algorithm. Currently no options are supported by the
base class.

=back

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Compute the options.
    my $topSize = $options{topSize} // 0;
    # Create the object.
    my $retVal = { topSize => $topSize };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Virtual Methods

=head3 type

    my $typeName = $compare->type();

Return the type string for this scoring method.

=cut

sub type {
    my ($self) = @_;
    return "bin" . $self->{topSize};
}


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
        push @scores, [$cv1->[$i], $cv2->[$i]];
    }
    my @sorted = sort { $b->[0] <=> $a->[0] } @scores;
    # The vector is now sorted by the CV1 scores, in order from largest to smallest.
    # We will now compare every two positions. If a subsequent position has a greater
    # CV2 score than a particular position, we have a failure and return 0.
    my $retVal = 1;
    $n = scalar(@sorted) - 1;
    if ($self->{topSize} && $self->{topSize} < $n) {
        $n = $self->{topSize};
    }
    # Loop through the pairs, searchign for a mismatch.
    for (my $i = 0; $i < $n && $retVal; $i++) {
        my $sortedI = $sorted[$i][1];
        my $sortedJ = $sorted[$i+1][1];
        if ($sortedJ > $sortedI) {
            $retVal = 0;
        }
    }
    # Return the indicator.
    return $retVal;
}


1;