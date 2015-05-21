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


package CPscore::Vector;

    use strict;
    use warnings;
    use base qw(CPscore);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs as a dot product
of similarity vectors. In this case, the similarity vector contains the best percent identity BLAST score for each
reference genome.

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
    # Get the identity score.
    my $iden = $sim->iden;
    # Update the score.
    my $retVal = ($iden > $oldScore ? $iden : $oldScore);
    # Return the result.
    return $retVal;
}

=head3 fix_vector

    my $keepFlag = $scoring->fix_vector($contig);

Adjust a similarity vector. This method performs post-processing on an assembled vector and allows the
scoring process one last chance to reject it. Note this method is called BEFORE the vector is normalized.

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
    my $v = $contig->vector;
    return sims_ok($v);
}



=head2 Utility Methods


=head3 sims_ok

Return TRUE if the incoming vector has a valid similarity configuration, else FALSE.

=cut

sub sims_ok
{
    my ($v) = @_;

    my $tot = 0;
    foreach $_ (@$v) { $tot += $_ }
    return (( $tot > 30) && ($tot < 10000))
}

1;