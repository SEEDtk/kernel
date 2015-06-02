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


package CPscore::Compare;

    use strict;
    use warnings;

=head1 Commmunity Pipeline Scoring Base Class

This is the base class for the community pipeline scoring objects that compute
vector comparison scores.

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
    # Create the object.
    my $retVal = {};
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
    return ref $self;
}

=head3 sortNeeded

    my $flag = $compare->sortNeeded();

Return TRUE if this is a normal scoring method that requires sorting, or FALSE if this is a binary
scoring method that returns only C<1> and C<0>. The default is TRUE.

=cut

sub sortNeeded {
    return 1;
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
    die "compare must be overridden.";
}


1;