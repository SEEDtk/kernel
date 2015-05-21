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


package CPscore;

    use strict;
    use warnings;
    use CPscore::Compare::Best;
    use CPscore::Compare::Binning;
    use CPscore::Compare::Ranking;
    use CPscore::Compare::Distance;
    use CPscore::Compare::DotProduct;

=head1 Community Pipeline Scoring

This is the base class for all community pipeline scoring objects. Its purpose is to compute the scoring
vector for each sample contig (one score per reference genome per contig) as well as the similarity score
for each contig pair.

This object has the following fields.

=over 4

=item compare

A L<Compare> object that can be used to compute the comparison score between two scoring vectors.

=item normalize

TRUE if the scoring vectors should be normalized, else FALSE.

=back

=head2 Special Methods

=head3 new

    my $scoring = CPscore->new($compareType, %options);

Create a new scoring object with the specified options.

=over 4

=item compareType

Comparison algorithm to use.

=item options

A hash containing the option values. The following keys are supported.

=over 8

=item normalize

TRUE if the scoring vectors should be normalized, else FALSE.

=back

=back

=cut

sub new {
    # Get the parameters.
    my ($class, $compareType, %options) = @_;
    # Get the comparison algorithm.
    my $compare;
    if ($compareType eq 'best') {
        $compare = CPscore::Compare::Best->new(%options);
    } elsif ($compareType eq 'bin') {
        $compare = CPscore::Compare::Binning->new(%options);
    } elsif ($compareType =~ /^bin(\d+)$/) {
        $compare = CPscore::Compare::Binning->new(topSize => $1);
    } elsif ($compareType eq 'dist') {
        $compare = CPscore::Compare::Distance->new(%options);
    } elsif ($compareType eq 'rank') {
        $compare = CPscore::Compare::Ranking->new(%options);
    } elsif ($compareType eq 'dot') {
        $compare = CPscore::Compare::DotProduct->new(%options);
    } else {
        die "Invalid comparison method $compareType.";
    }
    # Create the object.
    my $retVal = { compare => $compare, normalize => ($options{normalize} // 0) };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Utility Methods

=head3 type

    my $typeName = $scoring->type();

Return the type string for this scoring method.

=cut

sub type {
    my ($self) = @_;
    return join('-', ref($self), $self->{compare}->type, ($self->{normalize} ? 'n' : 'v'));
}

=head3 update_score

    $scoring->update_score($contig, $sim, $refID);

Process a similarity and use it to update the scoring vector coordinates.

=over 4

=item contig

A L<SampleContig> object containing the data for the contig.

=item sim

A L<Sim> object containing the results of a significant BLAST hit between the contig and a
reference genome.

=item refID

The ID of the target reference genome.

=back

=cut

sub update_score {
    # Get the parameters.
    my ($self, $contig, $sim, $refID) = @_;
    # Get the old score.
    my $oldScore = $contig->Hit($refID);
    # Compute the new score (this calls through to the subclass).
    my $score = $self->compute_score($contig, $sim, $refID, $oldScore);
    # Store it in the contig object.
    $contig->SetScore($refID, $score);
}

=head3 compare

    my $score = $scoring->compare($contig1, $contig2);

Compute the similarity score for a pair of contigs from their similarity vectors.

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
    # Call the comparison object.
    my $retVal = $self->{compare}->compare($contig1, $contig2);
    # Return the result.
    return $retVal;
}

=head3 adjust_vector

    my $keepFlag = $scoring->adjust_vector($contig);

Adjust a similarity vector. This method performs post-processing on an assembled vector and allows the
scoring process one last chance to reject it.

=over 4

=item contig

L<SampleContig> object whose vector is to be adjusted.

=item RETURN

Returns TRUE if the vector should be kept, FALSE if it should be discarded.

=back

=cut

sub adjust_vector {
    # Get the parameters.
    my ($self, $contig) = @_;
    # Call the subclass to fix the vector.
    my $retVal = $self->fix_vector($contig);
    # Normalize if necessary.
    if ($retVal && $self->{normalize}) {
        my $v = $contig->vector;
        my $norm = 0;
        # Note the vectors are sparse, so there is a real performance gain by skipping 0 positions.
        for my $x (@$v) {
            if ($x) {
                $norm += $x*$x;
            }
        }
        $norm = sqrt($norm);
        if ($norm > 0) {
            my $n = scalar @$v;
            for (my $i = 0; $i < $n; $i++) {
                $v->[$i] /= $norm;
            }
        }
    }
    # Return the keep/discard indicator.
    return $retVal;
}


=head2 Virtual Methods

=head3 compute_score

    my $score = $scoring->compute_score($contig, $sim, $refID, $oldScore);

Process a similarity and use it to update the scoring vector coordinates for the relevant genome.

=over 4

=item contig

A L<SampleContig> object containing the data for the contig.

=item sim

A L<Sim> object containing the results of a significant BLAST hit between the contig and a
reference genome.

=item refID

The ID of the target reference genome.

=item oldScore

The current score for the similarity between the contig and the reference genome.

=back

=cut

sub compute_score {
    # Get the parameters.
    my ($self, $contig, $sim, $refID, $oldScore) = @_;
    die "compute_score must be overridden.";
}

=head3 adjust_scores

    $scoring->adjust_scores($contig);

Perform a final adjustment on the scores for the contig before forming them into a vector. Use this method
for adjustments that are related solely to the scores themselves. For adjustments that are relevant to the
vector nature of the scores, use L</adjust_vector>. The default method does nothing.

=over 4

=item contig

L<SampleContig> object whose scores are to be adjusted.

=back

=cut

sub adjust_scores {
    # Get the parameters.
    my ($self, $contig) = @_;
}

=head3 fix_vector

    my $keepFlag = $scoring->fix_vector($contig);

Adjust a similarity vector. This method performs post-processing on an assembled vector and allows the
scoring process one last chance to reject it. The default method makes no adjustment and returns TRUE.

=over 4

=item contig

L<SampleContig> object whose vector is to be adjusted.

=item RETURN

Returns TRUE if the vector should be kept, FALSE if it should be discarded.

=back

=cut

sub fix_vector {
    # Get the parameters.
    my ($self, $contig) = @_;
    return 1;
}


1;