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


package SeedInfo;

    use strict;
    use warnings;
    use base qw(SEEDClient);


=head1 SEED Client Interface

This object provides an interface to the L<SEEDClient> object. Its main purpose is to abstract away
the need to know URLs. Simply use the L</new> method to create an object connected to the appropriate
SEED, and then call whichever SEEDClient method is required. For example, the following fragment would
get a hash of all the functional assignments for the pegs in a PubSEED genome.

    use SeedInfo;

    my $seed = SeedInfo->new("pubseed");
    my $fidH = $seed->get_genome_features([$genomeID], 'peg');
    my $funH = $seed->get_function($fidH->{$genomeID});

Because they use web services, all of the SEEDClient functions are batch oriented, and usually return
hash tables. The L<SEEDClient/get_genome_features> call returns a hash mapping the single incoming genome
ID to a list of feature IDs. The L<SEEDClient/get_function> call takes as input a list of feature IDs and
returns a hash mapping those feature IDs to functional assignments.

=head2 Special Methods

=head3 new

    my $seed = SeedInfo->new($type);

Return a L<SEEDClient> object for the named SEED instance.

=over 4

=item type

A string indicating the SEED instance. C<pubseed> for the PubSEED or C<core> for the Core SEED.

=cut

sub new {
    # Get the parameters.
    my ($class, $type) = @_;
    # Create the desired SEEDClient instance.
    my $retVal = SEEDClient::new($class, "http://$type.theseed.org/FIG/seed_svc");
    # Return it.
    return $retVal;
}

1;