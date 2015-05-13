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


package CPscore::Signal::Best;

    use strict;
    use warnings;
    use base qw(CPscore::Signal);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on a
conformity of the top-ranked signal strengths. The similarity score is I<1/(1+x)>, where I<x> is the percent
difference between the top coordinates in each vector. If the top coordinates are different, the
similarity score is 0.

In addition to the fields in the base class, this object contains the following fields.

=over 4

=item maxH

Reference to a hash mapping each contig ID to the position of its vector's highest coordinate.

=back

=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    my $retVal = 'signalbest';
    return $retVal;
}

=head3 init_vector_hash

    $scoring->init_vector_hash();

Initialize the scoring process. This process is used to compute the scoring vector
coordinates. For this object we initialize both hashes.

=cut

sub init_vector_hash {
    my ($self) = @_;
    $self->{vectorHash} = {};
    $self->{maxH} = {};
}


=head3 vector_compare

    my $score = $scoring->vector_compare($contig1, $cv1, $contig2, $cv2);

Compute the similarity score for a pair of contigs from their similarity vectors. Each similarity
vector contains the scores computed by L<CPscore::Signal/update_vector_hash> in order by the relevant reference
genome ID and formatted by L</store_vector>.

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
    my $retVal = 0;
    # Check the two vectors for a coordinate match.
    my $i1 = $self->{maxH}{$contig1};
    my $i2 = $self->{maxH}{$contig2};
    if ($i1 == $i2) {
        # The coordinates match, so compare them.
        my $v1 = $cv1->[$i1];
        my $v2 = $cv2->[$i2];
        if ($v1 > $v2) {
            $retVal = 1 / (1 + ($v1 - $v2)/$v1);
        } else {
            $retVal = 1 / (1 + ($v2 - $v1)/$v2);
        }
    }
    # Return the tally.
    return $retVal;
}

=head3 store_vector

    $scoring->store_vector($contig_vecs, $contig, $v);

Store a similarity vector in the vector hash.

=over 4

=item contig_vecs

Reference to a hash keyed on contig ID that is used to store the similarity vectors.

=item contig

ID of the contig whose vector is being stored.

=item v

Similarity vector to store.

=back

=cut

sub store_vector {
    # Get the parameters.
    my ($self, $contig_vecs, $contig, $v) = @_;
    # Determine whether or not the vector is acceptable. Also, save its max value coordinate.
    my $tot = 0;
    my $maxI;
    my $maxV = 0;
    my $n = scalar @$v;
    for (my $i = 0; $i < $n; $i++) {
        my $vi = $v->[$i];
        $tot += $vi;
        if ($vi > $maxV) {
            $maxI = $i;
            $maxV = $vi;
        }
    }
    if ($tot >= 500) {
        $contig_vecs->{$contig} = $v;
        $self->{maxH}{$contig} = $maxI;
    }
}



1;