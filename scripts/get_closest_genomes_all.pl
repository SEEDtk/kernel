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

=head1 Test the ability of get_closest_coreSEED_genomes.pl

    test_get_closest_genomes.pl > best.picks.from.coreSEED

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

my $min_hits = 100;  ### minimum number of kmer hits

my @tuples = $shrub->GetAll("Genome","",[],"Genome(id) Genome(name)");

foreach my $tuple (@tuples)
{
    my($g,$gs) = @$tuple;
    my $contig_file = &locate_contigs($shrub,$g);
    my @hits = map { (($_ =~ /^(\d+), (\S+)/) && ($1 >= $min_hits)) ? [$1,$2] : () }
               `perl get_closest_coreSEED_genomes.pl < $contig_file`;
    if (@hits == 0)
    {
	print STDERR join("\t",('could not be placed',$g,$gs)),"\n";
    }
    else
    {
	print join("\t",($hits[0]->[0],$hits[0]->[1],$g,$gs)),"\n";
	foreach my $hit (@hits)
	{
	    print join("\t",@$hit),"\n";
	}
	print "//\n";
    }
}

sub locate_contigs {
    my($schrub,$genome) = @_;

    return "a genome file";
}
    
