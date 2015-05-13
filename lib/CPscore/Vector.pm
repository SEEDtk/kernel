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
reference genome. The vector can be presented as is or normalized to a length of 100.

In addition to the fields in the base class, this object has the following fields.

=over 4

=item normalize

TRUE if the vectors should be normalized, else FALSE.

=back

=head2 Special Methods

=head3 new

    my $scoring = CPscore::Vector->new(%options);

Create a new scoring object with the specified options.

=over 4

=item options

A hash containing the option values. The following keys are supported.

=over 8

=item normalize

If TRUE, then each similarity vector will be normalized to a length of 1 before it is used to compute the
similarity score. The default is FALSE.

=back

=back

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Create the object.
    my $retVal = CPscore::new($class, %options);
    my $normalize = ($options{normalize} ? 1 : 0);
    $retVal->{normalize} = $normalize;
    # Return it.
    return $retVal;
}


=head2 Virtual Methods

=head3 type

    my $typeName = $scoring->type();

Return the type sttring for this scoring method.

=cut

sub type {
    my ($self) = @_;
    my $retVal = ($self->{normalize} ? 'normvec' : 'vector');
    return $retVal;
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
    # Get the identity score.
    my $iden = $sim->iden;
    # Get the vector hash.
    my $vecH = $self->{vectorHash};
    # Update the score.
    $vecH->{$contigID}{$refID} //= $iden;
    if ($iden > $vecH->{$contigID}{$refID}) {
        $vecH->{$contigID}{$refID} = $iden;
    }
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
    my ($self, $contig1, $cv1, $contig2, $cv2) = @_;
    # Return the dot product of the vectors.
    my $retVal = dot_product($cv1, $cv2);
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
    # Determine whether or not the vector is acceptable.
    if (sims_ok($v)) {
        # Normalize it if necessary.
        if ($self->{normalize}) {
            $v = unit_vector($v);
        }
        # Store it in the hash.
        $contig_vecs->{$contig} = $v;
    }
}


=head2 Utility Methods

=head3 dot_product

Compute the dot product of two vectors.

=cut

sub dot_product
{
    my ( $v1, $v2 ) = @_;

    my $tot = 0;
    my $i;
    for ( $i = 0 ; ( $i < @$v1 ) ; $i++ )
    {
        if ( $v1->[$i] && $v2->[$i] )
        {
            $tot += $v1->[$i] * $v2->[$i];
        }
    }
    return $tot;
}

=head3 unit_vector

Compute the normalized version of a vector.

=cut

sub unit_vector
{
    my ($v) = @_;

    my $tot = 0;
    my $uv  = [];
    my $i;
    for ( $i = 0 ; ( $i < @$v ) ; $i++ )
    {
        my $x = $v->[$i];
        if ( defined($x) )
        {
            $tot += $x * $x;
        }
    }

    my $nf = sqrt($tot);
    for ( $i = 0 ; ( $i < @$v ) ; $i++ )
    {
        my $x = $v->[$i];
        $x = $x ? $x : 0;
        if ($nf)
        {
            my $y = $x / $nf;
            if ( $y > 1 ) { $y = 1 }
            push( @$uv, sprintf( "%0.2f", $y ) );
        }
        else
        {
            push( @$uv, 0 );
        }
    }
    return $uv;
}

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