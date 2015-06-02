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


package SampleContig;

    use strict;
    use warnings;

=head1 Community Sample Contig Descriptor

This object represents a contig from a community sample. It contains the contig ID, its coverage amount, its
length, its universal role composition, and its match score for each reference genome.

The object contains the following fields.

=over 4

=item id

The contig ID.

=item len

The contig length, in base pairs.

=item covg

The contigs coverage value, an indication of how many times the contig was read by the sequencing process.

=item scores

Reference to a hash that maps each reference genome ID to a score that indicates its similarity to this contig.

=item vector

Reference to a list containing the contig's scoring vector. This is only filled in after a call to L</Score>.

=item bestCoord

Position in the scoring vector of the largest score. This is only filled in after a call to L</Score>.

=item roles

Reference to a hash listing the universal roles in this contig.

=item hits

Reference to a hash mapping each reference genome ID to the number of times it scored a hit on this contig.

=back

=head2 Special Methods

=head3 new

    my $contig = SampleContig->new($contigID, %options);

Create a new contig object.

=over 4

=item contigID

The ID of the contig.

=item options

A hash of options, containing zero or more of the following keys.

=over 8

=item len

The length of the contig, in base pairs. If this option is omitted, the length is computed from the contig ID.

=item covg

The coverage value of the contig. If this option is omitted, the coverage is computed from the contig ID.

=back

=back

=cut

sub new {
    # Get the parameters.
    my ($class, $contigID, %options) = @_;
    # Compute the length and coverage.
    my $len = $options{len};
    if (! defined $len) {
        if ($contigID =~ /length_(\d+)/) {
            $len = $1;
        } else {
            die "Contig ID $contigID does not specify a length and no \"len\" option specified.";
        }
    }
    my $covg = $options{covg};
    if (! defined $covg) {
        if ($contigID =~ /covg?_(\d+(?:\.\d+)?)/) {
            $covg = $1;
        } else {
            die "Contig ID $contigID does not specify a coverage and no \"covg\" option specified.";
        }
    }
    # Create the object.
    my $retVal = {  id =>       $contigID,
                    len =>      $len,
                    covg =>     $covg,
                    scores =>   {},
                    roles =>    {},
                    hits =>     {}
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Static Methods

=head3 get_contig

    my $contig = get_contig(\%contigs, $contigID);

Return the specified contig from the specified hash of L<SampleContig> objects. If the contig does not exist,
create a new object for it.

=over 4

=item contigs

Reference to a hash of L<SampleContig> objects, keyed by contig ID.

=item contigID

ID of the sample contig desired.

=item RETURN

Returns a L<SampleContig> object for the specified contig.

=back

=cut

sub get_contig {
    # Get the parameters.
    my ($contigs, $contigID) = @_;
    # Look for the contig.
    my $retVal = $contigs->{$contigID};
    # Does it exist?
    if (! defined $retVal) {
        # No, create it.
        $retVal = SampleContig->new($contigID);
        # Store it for future retrieval.
        $contigs->{$contigID} = $retVal;
    }
    # Return the object found.
    return $retVal;
}

=head2 Public Manipulation Methods

=head3 SetScore

    $contig->SetScore($genomeID, $score);

Store the similarity score for this contig relative to the specified reference genome.

=over 4

=item genomeID

ID of the relevant reference genome.

=item score

Value of the score to store.

=back

=cut

sub SetScore {
    # Get the parameters.
    my ($self, $genomeID, $score) = @_;
    # Store the score in the scoring hash.
    $self->{scores}{$genomeID} = $score;
}

=head3 SetRole

    $contig->SetRole($role);

Denote that the specified role is found in this contig.

=over 4

=item role

Name of the role found.

=back

=cut

sub SetRole {
    # Get the parameters.
    my ($self, $role) = @_;
    # Store the role in the role hash.
    $self->{roles}{$role} = 1;
}

=head3 Score

    my $okFlag = $contig->Score(\@genomes);

Form the scoring vector for this contig. The scoring vector will consist of the scores for the specified genomes,
in the order presented.

=over 4

=item genomes

Reference to a list of reference genome IDs.

=item RETURN

Return TRUE if at least one non-zero coordinate was found, else FALSE.

=back

=cut

sub Score {
    # Get the parameters.
    my ($self, $genomes) = @_;
    # Get the scoring hash.
    my $scores = $self->{scores};
    # Create the scoring vector.
    my @vector = map { $scores->{$_} // 0 } @$genomes;
    # Store it.
    $self->{vector} = \@vector;
    # Compute the best coordinate.
    my ($best, $bestVal) = (0, $vector[0]);
    for (my $i = 1; $i < scalar(@vector); $i++) {
        if ($vector[$i] > $bestVal) {
            $best = $i;
            $bestVal = $vector[$i];
        }
    }
    $self->{bestCoord} = $best;
    # Return TRUE if $bestVal is nonzero.
    return ($bestVal > 0);
}

=head3 Hit

    my $score = $contig->Hit($genomeID);

Record a hit by the specified reference genome against this contig and return the genome's current
similarity score.

=over 4

=item genomeID

ID of the relevant reference genome.

=item RETURN

Returns the current similarity score for the given reference genome.

=back

=cut

sub Hit {
    # Get the parameters.
    my ($self, $genomeID) = @_;
    # Record the hit.
    $self->{hits}{$genomeID}++;
    # Return the score.
    return $self->ref_score($genomeID);
}

=head2 Query Methods

=head3 ref_score

    my $score = $contig->ref_score($genomeID);

Return the score for the specified reference genome's similarity to this contig.

=over 4

=item genomeID

ID of the relevant reference genome.

=item RETURN

Returns the current score for the indicated reference genome.

=back

=cut

sub ref_score {
    my ($self, $genomeID) = @_;
    return $self->{scores}{$genomeID} // 0;
}

=head3 ref_hit

    my $hitCount = $contig->ref_hit($genomeID);

Return the number of hits against this contig by the specified reference genome.

=over 4

=item genomeID

ID of the relevant reference genome.

=item RETURN

Returns the hit count for this contig and the relevant reference genome.

=back

=cut

sub ref_hit {
    my ($self, $genomeID) = @_;
    return $self->{hits}{$genomeID} // 0;
}

=head3 scores

    my $scoreH = $contig->scores;

Return the hash of similarity scores for this contig.

=cut

sub scores {
    return $_[0]->{scores};
}


=head3 vector

    my $vector = $contig->vector($i);

Return the scoring vector or a coordinate in the scoring vector. This method returns C<undef> if L</Score>
has not yet been called.

=over 4

=item i

Index of the coordinate whose scores is desired, or C<undef> if the entire vector is desired.

=item RETURN

Returns the relevant vector coordinate, or the entire vector as a list reference.

=back

=cut

sub vector {
    my ($self, $i) = @_;
    my $retVal = $self->{vector};
    if ($retVal && defined $i) {
        $retVal = $retVal->[$i];
    }
    return $retVal;
}

=head3 id

    my $id = $contig->id;

Return the contig ID.

=cut

sub id {
    return $_[0]->{id};
}

=head3 len

    my $len = $contig->len;

Return the contig length.

=cut

sub len {
    return $_[0]->{len};
}

=head3 covg

    my $covg = $contig->covg;

Return the contig coverage.

=cut

sub covg {
    return $_[0]->{covg};
}

=head3 bestCoord

    my $bestCoord = $contig->bestCoord;

Return the index of the best coordinate in the scoring vector. This method returns C<undef> if
L</Score> has not yet been called.

=cut

sub bestCoord {
    return $_[0]->{bestCoord};
}

=head3 roles

    my $roleHash = $contig->roles

Return the hash of universal roles for this contig. The hash is keyed on role name, and has a key for each role found
in the contig.

=cut

sub roles {
    return $_[0]->{roles};
}


1;