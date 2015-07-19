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


package CPscore::Basis;

    use strict;
    use warnings;

=head1 Basis Vector Computation

This object is used to compute the basis vector for a community sample scoring operation. The basis vector
is the list of reference genomes whose coordinates are to be used. The default is to use all of them. Subclasses
can override the L</compute> method to propose an alternative algorithm.

=head2 Special Methods

=head3 new

    my $basis = CPscore::Basis->new(%options);

Construct a new basis-computation object.

=over 4

=item options

A hash of options. Currently there are none.

=back

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Create the object.
    my $retVal = {};
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Virtual Methods

=head3 type

    my $typeString = $basis->type;

Return a string description of this basis computation method.

=cut

sub type {
    return 'normal';
}

=head3 compute

    my $genomeList = $basis->compute(\%contigHash, \@genomeIDs);

Compute the list of reference genomes to use as the basis for the scoring vectors. This will be a subset
of the incoming genome ID list.

=over 4

=item contigHash

Reference to a hash of L<SampleContig> objects containing the scoring data for all the sample contigs,
keyed by contig ID.

=item genomeIDs

Reference to a list of genome IDs. Only genome IDs in this list can be part of the basis.

=item RETURN

Returns a list of genome IDs whose scores should be represented in the scoring vectors.

=back

=cut

sub compute {
    # Get the parameters.
    my ($self, $contigHash, $genomeIDs) = @_;
    # Return the incoming list of genome IDs;
    my @retVal = @$genomeIDs;
    print STDERR scalar(@retVal) . " genomes returned.\n";
    return \@retVal;
}

1;