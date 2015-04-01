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
use File::Copy::Recursive;

=head1 Project Subsystems Onto Genomes

    solid_projection_pipeline [options] --subsystem=<sub_id> --dataD=<output directory> < genome_list

Look at the features used by a subsystem in a particular set of genomes, and use the information
to project the subsystem onto as many other genomes as possible.

The program must first compute a set of parameters for determining valid instances of a role. These
parameters are derived by looking at genomes which have solid
instances of the subsystem in question in the core SEED. A solid instance generally has a complete
set of roles for its variant and only one peg implementing each role. For a given role, we compare
the lengths of the implementing protein sequences to determine a mean and standard deviation. In
addition, we BLAST the proteins against each other to compute a similarity threshhold we can use
to decide if a protein is a valid instance of a function. This gives us a set of criteria for
projecting the role onto a new genome.

Once the parameters are known, we look at other genomes in the database and look at pegs whose
functions contain the roles in our target subsystem. If the peg's protein meets the criteria
determined in the first phase, then we consider it a match for the subsystem role. If a genome
has enough role matches to fill a variant, we recommend a projection of the subsystem
onto it.

=head2 Output

The output from this method is stored in the data directory specified in the command-line
options. The files C<parms.1> and C<parms.2> contain the computed parameters for the
projection in JSON format. The file C<projections.1> contains the core SEED genomes predicted
as likely candidates to contain the subsystem as a result of the first pass over the data,
using the parameters computed from the incoming genomes.  These new genomes are combined with
the incoming genomes to compute a new set of parameters (this is what goes into C<parms.2>).
The final projections using these new parameters are written to C<projections.2>.

Each projections file is divided into sections terminated by a C<//> line. Each section represents
a prediction that a subsystem has an implementation in a genome. The first record of the section
has four tab-delimited columns as follows.

=over 4

=item 1

ID of the subsystem.

=item 2

ID of the genome that is predicted to hold an implementation of it.

=item 3

The variant code we believe the genome uses. If the genome does not appear to implement
the subsystem, this column will read C<not-active>.

=item 4

The ID of a template genome containing a nearly identical set of roles, that is associated
with the variant code.

=back

Next, if we have an active variant,  there is a series of one or more records describing
the pegs that implement roles in the subsystem. The series is terminated by a line of four
hyphens (C<---->). Each record in the series has three tab-delimited columns as follows.

=over 4

=item 1

ID of the peg that is believed to implement the role.

=item 2

ID of the role.

=item 3

ID of the function assigned to the peg.

=back

Finally, there is a series of records representing pegs that are assigned to functions
implementing roles for the subsystem, but which violated one or more of the parameters
determined by the analysis of the solid instances. Each record in this series is
multiple tab-delimited columns, the first being the peg ID, and each additional column
being text that describes why the peg is considered problematic.

=head2 Parameters

There are no positional parameters. The standard input is presumed to contain a list of
genomes that are considered solid examples of the subsystem. The command-line options are
those found in L<Shrub/script_options> and L<ScriptUtils/ih_options> plus the following.

=over 4

=item subsystem

The ID of a subsystem to be projected.

=item dataD

Data directory used to store the match parameters and the results.

=item privilege

The required privilege level for annotations to be considered when projecting.

=item chooseGenomes

If this option is specified, the initial list of genomes will be computed automatically
rather than being read from the input.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'subsystem|s=s', 'ID of the subsystem to process', { required => 1 } ],
    [
        'dataD|d=s',
        'Data Directory for Subsystem Projection',
        { required => 1 }
    ],
    ['chooseGenomes|G', 'choose the input genomes automatically'],
    [ 'privilege=i', 'privilege level of roles', { default => Shrub::PRIV } ]
);
my $subsystem_id = $opt->subsystem;
my $dataD        = $opt->datad;
my $privilege    = $opt->privilege;
my $shrub        = Shrub->new_for_script($opt);

# Insure we have a data directory.
if (! -d $dataD) {
    File::Copy::Recursive::pathmk($dataD);
}

my @genomes;
# Are we choose the genomes automatically?
if ($opt->choosegenomes) {
    # Yes. Search the database for appropriate genomes.
    @genomes = Projection::choose_genomes($shrub, $subsystem_id);
} else {
    my $ih       = ScriptUtils::IH( $opt->input );
    @genomes  = map { ( $_ =~ /(\d{3,10}\.\d+)/ ) ? $1 : () } <$ih>;
}
my @full_set = @genomes;

# Compute the parameters for the projection.
my $parms =
  Projection::compute_properties_of_solid_roles( $shrub, $subsystem_id,
    \@genomes );

# Save them to disk.
SeedUtils::write_encoded_object( $parms, "$dataD/parms.1" );

# Get the remaining core genomes in the database.
my %small_set = map { ( $_ => 1 ) } @genomes;
my @to_call =
  grep { !$small_set{$_} }
  $shrub->GetFlat( 'Genome', 'Genome(core) = ?', [1], 'id' );

# Project the subsystem onto them.
open( my $oh, ">$dataD/projections.1" )
  || die "Could not open first projection file: $!";
my @found =
  Projection::project_solid_roles( $shrub, $subsystem_id, $privilege, \@to_call,
    $parms, $oh );
close $oh;
undef $oh;

# Add the new genomes to our set to process.
push @full_set, @found;

# Compute a new set of parameters.
my $parms2 =
  Projection::compute_properties_of_solid_roles( $shrub, $subsystem_id,
    \@full_set );

# Save the new parms to disk.
SeedUtils::write_encoded_object( $parms2, "$dataD/parms.2" );

# Get a list of all the genomes.
my @all = grep { !$small_set{$_} } $shrub->GetFlat( 'Genome', '', [], 'id' );

# Project the subsystem with the new parameters.
open( $oh, ">$dataD/projections.2" )
  || die "Could not open second projection file: $!";
Projection::project_solid_roles( $shrub, $subsystem_id, $privilege, \@to_call,
    $parms2, $oh );

