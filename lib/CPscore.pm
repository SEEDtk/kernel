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

=head1 Community Pipeline Scoring

This is the base class for all community pipeline scoring objects. Its purpose is to compute the scoring
vector for each sample contig (one score per reference genome per contig) as well as the similarity score
for each contig pair.

This object has the following fields.

=over 4

=item vectorHash

A reference to a two-dimensional hash, keyed on sample contig ID followed by reference genome ID. The hash tracks
the similarity score for each sample contig / reference genome pair.

=back

=head2 Special Methods

=head3 new

    my $scoring = CPscore->new(%options);

Create a new scoring object with the specified options.

=over 4

=item options

A hash containing the option values. No keys are supported for the base class.

=back

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Create the object.
    my $retVal = { vectorHash => {} };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    die "type method not overridden.";
}

=head3 init_vector_hash

    $scoring->init_vector_hash();

Initialize the scoring process. This process is used to compute the scoring vector
coordinates. For this object we simply create an empty hash.

=cut

sub init_vector_hash {
    my ($self) = @_;
    $self->{vectorHash} = {};
}

=head3 get_vector_hash

    $scoring->get_vector_hash()

Return the hash containing the scoring vector coordinates.

=cut

sub get_vector_hash {
    my ($self) = @_;
    return $self->{vectorHash};
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
    die "update_vector_hash must be overridden.";
}

=head3 vector_compare

    my $score = $scoring->vector_compare($contig1, $cv1, $contig2, $cv2);

Compute the similarity score for a pair of contigs from their similarity vectors. Each similarity
vector contains the scores computed by L</update_vector_hash> in order by the relevant reference
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
    my ($self, $contig1, $cv1, $contig2, $cv2) = @_;    die "vector_compare must be overridden.";
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
    # Determine whether or not the vector is acceptable.
    my $tot = 0;
    for my $vi (@$v) {
        $tot += $vi;
    }
    if ($tot >= 500) {
        $contig_vecs->{$contig} = $v;
    }
}




1;