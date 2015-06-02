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


package GenomeFinder;

    use strict;
    use warnings;
    use Shrub::GTO;
    use SeedURLs;
    use Shrub;
    use LWP::Simple;
    use GenomeTypeObject;
    use JSON::XS;

=head1 Create a GenomeTypeObject

This library contains methods that will create a L<GenomeTypeObject> for a genome from one of several SEED sources.
It will first attempt to create the GTO from the Shrub database. If this fails, it will query the various online
SEED servers.

=head2 Public Methods

=head3 find

    my $gto = GenomeFinder::find($genomeID, $shrub);

Return the genome with the specified ID.

=over 4

=item genomeID

ID of the desired genome.

=item shrub

L<Shrub> database object to search for the genome. If omitted, one will be created.

=item RETURN

Returns a L<GenomeTypeObject> for the specified genome, if one exists, and C<undef> otherwise.

=back

=cut

# List of the seed server types to check for the genome, in order from last to first.
use constant SEED_TYPES => ['pseed', 'pubseed'];

sub find {
    # Get the parameters.
    my ($genomeID, $shrub) = @_;
    # Insure we have a Shrub database.
    $shrub //= Shrub->new();
    # Check it for the genome.
    my $retVal = Shrub::GTO->new($shrub, $genomeID);
    # Loop through the seed servers until we find it.
    my @servers = @{SEED_TYPES()};
    while (! $retVal && @servers) {
        my $server = pop @servers;
        $retVal = seed($server, $genomeID);
    }
    # Return the GTO found.
    return $retVal;
}

=head3 seed

    my $gto = GenomeFinder::seed($serverName, $genomeID);

Locate a genome on a specific SEED server.

=over 4

=item serverName

The name (or URL) of the SEED server containing the genome.

=item genomeID

The ID of the desired genome.

=item RETURN

Returns a L<GenomeTypeObject> for the genome, or undefined if the genome was not found on the server.

=back

=cut

sub seed {
    # Get the parameters.
    my ($serverName, $genomeID) = @_;
    # Compute the URL.
    my $serverURL = SeedURLs::url($serverName);
    # Declare the return variable.
    my $retVal;
    # Ask for the GTO.
    my $json = LWP::Simple::get("$serverURL/genome_object.cgi?genome=$genomeID");
    # If we found it, decode it.
    if ($json) {
        $retVal = decode_json($json);
        GenomeTypeObject->initialize($retVal);
    }
    # Return the result.
    return $retVal;
}

1;