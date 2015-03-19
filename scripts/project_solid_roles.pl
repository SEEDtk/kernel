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
use FIG_Config;
use Shrub;
use ScriptUtils;
use Projection;
use Data::Dumper;

=head1 project solid instances of subsystems to new genomes

    project_solid_roles -s Subsystem-id < GenomesFile [ options ]

project to solid instances of a subsystem in new genomes.

=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -s Subsystem-ID

Note that this tool takes the short, unique id as input -- not the
full subsystem name.

=item -d encodedF

Names the file that contains the parameters for the projections of the subsystem.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s', 'ID of the subsystem to process' ],
    [ 'privilege|p=i', 'Privilege level of function and role names', { default => Shrub::PRIV } ],
    [ 'encodedF|e=s',     'parameter file for Subsystem Projection' ]
);
my $subsystem_id = $opt->subsystem;
my $eFile        = $opt->encodedf;
my $privilege    = $opt->privilege;

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

# Open the input file.
my $ih = ScriptUtils::IH( $opt->input );
my @genomes = map { ( $_ =~ /(\d+\.\d+)/ ) ? $1 : () } <$ih>;
close($ih);

my $parms = &Projection::read_encoded_object($eFile);

Projection::project_solid_roles($shrub, $subsystem_id, $privilege, \@genomes, $parms, \*STDOUT);
