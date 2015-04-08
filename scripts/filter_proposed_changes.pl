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
use Shrub;
use ScriptUtils;
use Data::Dumper;
use File::Copy::Recursive;

=head1 Construct List of Proposed additions to subsystems

    filter_proposed_changes -c Current-memberships -p proposed > condensed-changes

We have a pipeline for proposing projections of subsystems.  The report needs to be condensed
into a set of proposed updates to subsystems.  We take as input
two pieces of data: the report giving projections and the current memberships in subsystems.

We need to be careful to filter out proposals in which the current membership/vc already 
exists.
=head2 Output

The proposed changes from the projection report that have not already been made.

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

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    ScriptUtils::ih_options(),
    [ 'current|c=s', 'A table of current subsystem membersjips', { required => 1 } ],
    [ 'projections|p=s', 'A pipeline report of proposed projections', { required => 1 } ]
);
my $currentF     = $opt->current;
my $proposed      = $opt->projections;
my $shrub        = Shrub->new_for_script($opt);
