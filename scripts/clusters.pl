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
use Shrub;
use ScriptUtils;
use Data::Dumper;
use File::Slurp;
use SeedUtils;
use Clusters;

=head1 Compute Clusters

    clusters [-m MaxGap] [ options ] > Clusters

computes clusters on the chromosome of pegs from a subsystem (perhaps
multiple subsystems).  The primary output is the input table with 2 columns
appended:

    RowId   is the id of a row in a populated subsystem

    Cluster is a comma-separated set of peg ids
    
MaxGap is optionl.

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

-m MaxGapBetweenGenesInCluster

    the maximum gap between members of a cluster.  Defaults to 2000bp.

-c N

    the column containing genome IDs

=cut

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '',
                     Shrub::script_options(), ScriptUtils::ih_options(),
		     [ 'maxgap|m=i', { default => 2000 } ]
    );

my $max_gap    = $opt->maxgap;
my $shrub      = Shrub->new_for_script($opt);
####  Fails to pick up default   !!!!!
#     die $max_gap;
while (defined($_ = <STDIN>))
{
    chop;
    my $genome = $_;
    
    my $peg_to_subsys_row = &Clusters::subsys_mems($genome,$shrub);
    my @fidPairs   = grep { $_->[1] } map { [$_,$shrub->loc_of($_)] } keys(%$peg_to_subsys_row);
    my @sorted     = sort { BasicLocation::Cmp($a->[1], $b->[1]) } @fidPairs;
    my @clusters = &Clusters::cluster_subsys_fids($max_gap,\@sorted,$peg_to_subsys_row);
    foreach my $cluster (@clusters)
    {
	&print_cluster($cluster);
    }
}

sub print_cluster {
    my($cluster) = @_;

    my($row,$pegs) = @$cluster;
    print $row,"\t",join(",",@$pegs),"\n";
}


