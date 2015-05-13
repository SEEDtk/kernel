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


package CPscore::Signal;

    use strict;
    use warnings;
    use base qw(CPscore);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on
similarity signal vectors. In this case, the similarity vector contains the sum of the percent identity times the length
of all the BLAST matches to each reference genome. The resulting number is an indication of how much of the sample contig
matches each reference genome. This is a base class for several scoring methods that use signal strength vectors.

=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    die "type method must be overridden.";
}

=head3 update_vector_hash

    $scoring->update_vector_hash($sim, $refID);

Process a similarity and use it to update the scoring vector coordinates.

=over 4

=item sim

A L<Sim> object containing the results of a significant BLAST hit between a contig and a
reference genome.

=item refID

The ID of the target reference genome.

=back

=cut

sub update_vector_hash {
    # Get the parameters.
    my ($self, $sim, $refID) = @_;
    # Extract the contig ID.
    my $contigID = $sim->id2();
    # Get the match strength.
    my $strength = $sim->iden * abs($sim->e2() - $sim->b2()) / 100;
    # Get the vector hash.
    my $vecH = $self->{vectorHash};
    # Update the score.
    $vecH->{$contigID}{$refID} //= 0;
    $vecH->{$contigID}{$refID} += $strength;
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
    die "vector_compare must be overridden."
}

1;