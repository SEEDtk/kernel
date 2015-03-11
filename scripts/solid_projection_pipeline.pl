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
use Data::Dumper;

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '', Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s','ID of the subsystem to process'],
    [ 'dataD|d=s','Data Directory for Subsystem Projection']
);
my $subsystem_id = $opt->subsystem;
my $dataD        = $opt->datad;
my $shrub        = Shrub->new_for_script($opt);
$subsystem_id || die "missing subsystem id";
$dataD        || die "missing Data directory";

my $ih = ScriptUtils::IH($opt->input);
my @genomes = map { ($_ =~ /(\d{3,10}\.\d+)/) ? $1 : () } <$ih>;
my @full_set = @genomes;

open(GETPROP,"| perl properties_of_solid_roles.pl -s $subsystem_id -d $dataD/solid.projection.parms.1")
    || die "failed to run initial properies_of_solid_roles";
foreach my $g (@genomes)
{
    print GETPROP $g,"\n";
}
close(GETPROP);

my %small_set = map { ($_ => 1) } @genomes;
my @to_call = grep { ! $small_set{$_} } &get_genomes_from_database($shrub);


open(PROJ,"| perl project_solid_roles.pl -s $subsystem_id -d $dataD/projections.1")
    || die "failed projection 1";
foreach my $g (@to_call)
{
    print PROJ $g,"\n";
}
close(PROJ);

$/ = "\n//\n";
foreach my $called (`cat $dataD/projections.1`)
{
    if ($called =~ /^\S+\t(\S+)/)
    {
    push(@full_set,$1);
    }
}
$/ = "\n";

open(GETPROP,"| perl properties_of_solid_roles.pl -s $subsystem_id -d $dataD/solid.projection.parms.2")
    || die "failed to run second properies_of_solid_roles";
foreach my $g (@full_set)
{
    print GETPROP $g,"\n";
}
close(GETPROP);

open(PROJ,"| perl project_solid_roles.pl -s $subsystem_id -d $dataD/projections.2")
    || die "failed projection 2";
foreach my $g (@to_call)
{
    print PROJ $g,"\n";
}
close(PROJ);

sub get_genomes_from_database {
    my($shrub) = @_;

    my @tuples = $shrub->GetAll('Genome','',[],"Genome(id)");
    return map { $_->[0] } @tuples;
}

