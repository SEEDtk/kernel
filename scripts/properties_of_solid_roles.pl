#!/usr/bin/env perl
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

use strict;
use JSON::XS;
use warnings;
use Shrub;
use ScriptUtils;
use Projection;
use Data::Dumper;

=head1 compute properties of solid roles to be used in projections

    properties_of_solid_roles -s Subsystem-id < GenomesFile [ options ]

compute the properties of roles in a Subsystem using

=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -s Subsystem-ID

Note that this tool takes the short, unique id as input -- not the
full subsystem name.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '', Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s','ID of the subsystem to process', { required => 1}],
);

my $subsystem_id = $opt->subsystem;

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH( $opt->input );
# Read the genomes.
my @genomes = map { ( $_ =~ /(\d+\.\d+)/ ) ? $1 : () } <$ih>;
close($ih);
# Compute the properties.
my $parms = Projection::compute_properties_of_solid_roles($shrub, $subsystem_id, \@genomes);
# Write them to the standard output.
SeedUtils::write_encoded_object($parms, \*STDOUT);
