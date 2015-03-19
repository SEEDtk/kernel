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
use warnings;
use JSON::XS;
use FIG_Config;
use Shrub;
use ScriptUtils;
use Data::Dumper;
use Projection;

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '', Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s','ID of the subsystem to process', { required => 1 }],
    [ 'dataD|d=s','Data Directory for Subsystem Projection', { required => 1 }],
    [ 'privilege=i', 'privilege level of roles', { default => Shrub::PRIV }]
);
my $subsystem_id = $opt->subsystem;
my $dataD        = $opt->datad;
my $privilege    = $opt->privilege;
my $shrub        = Shrub->new_for_script($opt);

my $ih = ScriptUtils::IH($opt->input);
my @genomes = map { ($_ =~ /(\d{3,10}\.\d+)/) ? $1 : () } <$ih>;
my @full_set = @genomes;

# Compute the parameters for the projection.
my $parms = Projection::compute_properties_of_solid_roles($shrub, $subsystem_id, \@genomes);
# Save them to disk.
Projection::write_encoded_object($parms, "$dataD/parms.1");
# Get the remaining genomes in the database.
my %small_set = map { ($_ => 1) } @genomes;
my @to_call = grep { ! $small_set{$_} } &get_genomes_from_database($shrub);
# Project the subsystem onto them.
open(my $oh, ">$dataD/projections.1") || die "Could not open first projection file: $!";
my @found = Projection::project_solid_roles($shrub, $subsystem_id, $privilege, \@to_call, $parms, $oh);
close $oh; undef $oh;
# Add the new genomes to our set to process.
push @full_set, @found;
# Compute a new set of parameters.
my $parms2 = Projection::compute_properties_of_solid_roles($shrub, $subsystem_id, \@full_set);
# Save the new parms to disk.
Projection::write_encoded_object($parms2, "$dataD/parms.2");
# Project the subsystem with the new parameters.
open($oh, ">$dataD/projections.2") || die "Could not open second projection file: $!";
Projection::project_solid_roles($shrub, $subsystem_id, $privilege, \@to_call, $parms2, $oh);


sub get_genomes_from_database {
    my($shrub) = @_;

    my @tuples = $shrub->GetAll('Genome','',[],"Genome(id)");
    return map { $_->[0] } @tuples;
}

