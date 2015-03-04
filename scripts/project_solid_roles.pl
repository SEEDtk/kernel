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

=item -d dataD

Names the directory that retains state for the projections of the subsystem.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '', Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s','ID of the subsystem to process'],
    [ 'privilege|p=i','Privilege level of function and role names'],
    [ 'dataD|d=s','Data Directory for Subsystem Projection']
);
my $subsystem_id = $opt->subsystem;
my $dataD        = $opt->datad;
my $privilege    = $opt->privilege;
if (! $privilege) { $privilege = 2 }

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

# Open the input file.
my $ih = ScriptUtils::IH( $opt->input );
my @genomes = map { ( $_ =~ /(\d+\.\d+)/ ) ? $1 : () } <$ih>;
close($ih);

my $parms = &Projection::read_encoded_object("$dataD/solid.projection.parameters");

my @tuples = $shrub->GetAll("Subsystem2Role Role Role2Function Function Function2Feature Feature Feature2Contig",
                         "(Subsystem2Role(from-link) = ?) AND (Function2Feature(security) = $privilege)",
                         [$subsystem_id,$privilege],
             "Subsystem2role(to-link) Function2Feature(to-link) Function(description)
                          Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(dir)");
my %relevant;
foreach $_ (@tuples)
{
    my($role,$peg,$func,$contig,$beg,$strand) = @$_;
    my $g = &SeedUtils::genome_of($peg);
    $relevant{$g}->{$peg} = [$role,$func,[$contig,$beg,$strand]];
}

my $state = { (relevant => \%relevant) };

foreach my $g (@genomes)
{
    my $projection = &Projection::project_subsys_to_genome($shrub,$g,$subsystem_id,$state,$parms);
    if (my $vc = $projection->{vc})
    {
    print join("\t",($subsystem_id,$g,$vc)),"\n";
    my $calls = $projection->{calls};
        foreach my $call (@$calls)
    {
        my($peg,$role,$func) = @$call;
        print "\t",join("\t",($peg,$role,$func)),"\n";
    }
    print "//\n";
    }
}
