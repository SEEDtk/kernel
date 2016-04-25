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
use Shrub;
use ScriptUtils;

=head1 Display Database Statistics

    dbstats.pl [ options ] 

This script displays the number of genomes, features, taxonomic groups, roles, and subsystems in the
Shrub database. The output is tab-delimited, with an object type in the first column and a count in the
second. 

=head2 Parameters

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', 
        Shrub::script_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Loop through the object types.
my @types = qw(TaxonomicGrouping Genome Feature Role Subsystem Protein);
for my $type (@types) {
    my $count = $shrub->GetCount($type, '', []);
    print "$type\t$count\n";
}
