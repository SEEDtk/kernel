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
use Data::Dumper;

=head1 Get Roles that are universal and occur as singletons in prokaryotic genomes

    get_univ_singletons [ options ] > RolesAndCounts

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=cut

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '', Shrub::script_options(), ScriptUtils::ih_options(),
    [] );
my $ih = ScriptUtils::IH( $opt->input );

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my @tuples = $shrub->GetAll("Function Function2Feature Feature",
			    "Function2Feature(security) = ?",[2],
			    "Function(description) Function2Feature(to-link)");

my %funcH;
foreach my $tuple (@tuples)
{
    my($func,$fid) = @$tuple;
    if (! &SeedUtils::hypo($func))
    {
	if ($fid =~ /^fig\|(\d+\.\d+)\.peg\./)
	{
	    $funcH{$func}->{$1}++;
	}
    }
}
my %counts;
foreach my $func (keys(%funcH))
{
    my $genomesH = $funcH{$func};
    my @genomes = keys(%$genomesH);
    foreach my $g (keys(%$genomesH))
    {
	if ($genomesH->{$g} == 1) { $counts{$func}++ }
    }
}
	
my @best = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
foreach $_ (@best)
{
    print join("\t",($counts{$_},$_)),"\n";
}
