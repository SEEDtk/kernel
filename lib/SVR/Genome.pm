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


package SVR::Genome;

    use strict;
    use warnings;
    use Shrub::GTO;

=head1 Genome Server Module

This module takes as input a genome ID and returns a L<GenomeTypeObject> for that genome.

=head2 Special Methods

=head3 exec

    my $gto_json = SVR::Genome::exec($cgi, $shrub, \%parms);

Return a L<GenomeTypeObject> json string for a genome.

=over 4

=item cgi

The L<CGI> object for the current web request.

=item shrub

A L<Shrub> object for accessing the database.

=item parms

A reference to a hash containing the following parameters.

=over 8

=item id

A 1-tuple containing the ID of the desired genome.

=back

=item RETURN

Returns a json string representing the genome object.

=back

=cut

sub exec {
    my ($cgi, $shrub, $parms) = @_;
    # This will be the return value.
    my $retVal;
    # Get the genome ID.
    my $gid = $parms->{id};
    if (! $gid) {
        die "No genome ID specified.";
    } else {
        # Get the GTO.
        my $gto = Shrub::GTO->new($shrub, $gid->[0]);
        if (! $gto) {
            die "Genome $gid not found.";
        } else {
            # Converty it to a string.
            $gto->destroy_to_file(\$retVal);
        }
    }
    return $retVal;
}

1;