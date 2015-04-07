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
use FIG_Config;
use ScriptUtils;
use Shrub;
use Projection;
use File::Copy::Recursive;


=head1 Project Multiple Subsystems Onto Shrub Genomes

    solid_projection_engine [ options ] outputDirectory

This script looks at the solid implementations of subsystems in the Shrub database and uses the
information to project the subsystems onto the remaining genomes. It is expected it will run a
long time.

For each subsystem, we first find the genomes that contain solid implementations of the subsystem.
A genome contains a solid implementation if it has an active variant (B<needs-curation> is FALSE)
and each cell in the subsystem row contains only one peg.For a given role, we compare the lengths
of the implementing protein sequences to determine a mean and standard deviation. In addition,
we BLAST the proteins against each other to compute a similarity threshhold we can use to decide
if a protein is a valid instance of a function. This gives us a set of criteria for projecting
the role onto a new genome.

Once the parameters are known, we look at other genomes in the database and look at pegs whose
functions contain the roles in our target subsystem. If the peg's protein meets the criteria
determined in the first phase, then we consider it a match for the subsystem role. If a genome
has enough role matches to fill a variant, we recommend a projection of the subsystem onto it.

We perform this two-phase process twice-- once against core genomes to compute a larger set
of genomes with solid instances, and again (with this larger set as a base) against the
non-core genomes.

=head2 Output

The positional parameter for this command is a directory name. For each subsystem processed,
a subdirectory having the same name as the subsystem ID is created under this output
directory. Five files are stored in this directory.

=over 4

=item genomes.tbl

A tab-delimited file containing a list of the solid genomes for the subsystem.

=item solid.parms.json

A json file containing the parameters computed from the solid genomes.

=item core.parms.json

A json file containing the parameters computed from the solid genomes and the projected core genomes.

=item solid.proj.tbl

A tab-delimited file containing the projections from the solid genomes onto the core genomes.

=item core.proj.tbl

A tab-delimited file containing the projections from the core genomes onto the non-core genomes.

=back

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

The single positional parameter is the name of the output directory. A directory for each
subsystem will be created (or reused) under this directory.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item subsystem

The ID of a subsystem to process. This parameter can be specified multiple times to specify
multiple subsystems. Specify C<core> to process all core subsystems, C<all> to process all
subsystems in the database. The default is C<core>.

=item missing

If specified, only subsystems that do not already have an output directory will be processed.
The default is to replace the output files in any existing output directory.

=item privilege

Privilege level to use for the annotations on the non-core genomes. The default is C<0> (public).

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outputDirectory', Shrub::script_options(),
        ['subsystem|s=s@', 'ID of subsystem to process, "all" for all, or "core" for all core subsystems', { default => ['core'] }],
        ['missing|m', 'process only new subsystems'],
        ['privilege|p=i', 'priviliege level for second-pass annotations', { default => Shrub::PUBLIC }],
        );
# Verify the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "Output directory is required.";
} elsif (-f $outDir) {
    die "Output directory $outDir is invalid.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir);
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the list of subsystems to process.
my %subs;
for my $subSpec (@{$opt->subsystem}) {
    if ($subSpec eq 'all') {
        map { $subs{$_} = 1 } $shrub->GetFlat('Subsystem', '', [], 'id');
    } elsif ($subSpec eq 'core') {
        map { $subs{$_} = 1 } $shrub->GetFlat('Subsystem', 'Subsystem(privileged) = ?', [1], 'id');
    } else {
        $subs{$subSpec} = 1;
    }
}
# Get a list of all the genomes, mapping each genome ID to its core flag. We will use this hash to
# determine the core and non-core genomes.
my %genomes = map { $_->[0] => [$_->[1], $_->[2]] } $shrub->GetAll('Genome', '', [], 'id core name');
print scalar(keys %genomes) . " in database.\n";
# Loop through the subsystems.
my $totSubs = scalar(keys %subs);
my $subsRun = 0;
for my $sub (sort keys %subs) {
    $subsRun++;
    print "Processing $sub ($subsRun of $totSubs).\n";
    # This will be set to FALSE if we want to skip this subsystem.
    my $okSub = 1;
    # This will hold the solid genomes.
    my %solids;
    # Do we have an output directory for this subsystem?
    my $dataD = "$outDir/$sub";
    if (-d $dataD && $opt->missing) {
        print "Subsystem $sub already has data. Skipped.\n";
        $okSub = 0;
    } else {
        # Insure we have an output directory.
        if (! -d $dataD) {
            mkdir $dataD;
            print "Creating subsystem output directory $dataD.\n";
        }
        # Get the solid genomes.
        %solids = map { $_ => 1 } Projection::choose_genomes($shrub, $sub);
        if (! keys %solids) {
            print STDERR "No solid genomes found for $sub: subsystem skipped.\n";
            $okSub = 0;
        } else {
            print scalar(keys %solids) . " solid genomes found for $sub.\n";
            open(my $oh, ">$dataD/genomes.tbl") || die "Could not open genome output file for $sub: $!";
            for my $solid (sort keys %solids) {
                print $oh join("\t", $solid, $genomes{$solid}[1]) . "\n";
            }
        }
    }
    # Only proceed if our subsystem is OK to process.
    if ($okSub) {
        my @full_set = keys %solids;
        print "Computing solid parameters.\n";
        # Compute the parameters for the projection.
        my $parms =
          Projection::compute_properties_of_solid_roles( $shrub, $sub,
            [keys %solids] );
        # Save them to disk.
        SeedUtils::write_encoded_object( $parms, "$dataD/solid.parms.json" );
        # Get the remaining core genomes in the database.
        my @to_call =
          grep { !$solids{$_} && $genomes{$_}[0] } keys %genomes;
        # Project the subsystem onto them.
        open( my $oh, ">$dataD/solid.proj.tbl" )
          || die "Could not open solid projection file: $!";
        my @found =
          Projection::project_solid_roles( $shrub, $sub, Shrub::PRIV, \@to_call,
            $parms, $oh );
        close $oh;
        undef $oh;

        # Add the new genomes to our set to process.
        push @full_set, @found;
        print "Computing core parameters.\n";
        # Compute a new set of parameters.
        my $parms2 =
          Projection::compute_properties_of_solid_roles( $shrub, $sub,
            \@full_set );
        # Save the new parms to disk.
        SeedUtils::write_encoded_object( $parms2, "$dataD/core.parms.json" );
        # Get a list of all the genomes.
        my @all = grep { !$solids{$_} && !$genomes{$_}[0] } keys %genomes;
        # Project the subsystem with the new parameters.
        open( $oh, ">$dataD/core.proj.tbl" )
          || die "Could not open second projection file: $!";
        Projection::project_solid_roles( $shrub, $sub, $opt->privilege, \@to_call,
            $parms2, $oh );
        close $oh;
        undef $oh;
    }
}
