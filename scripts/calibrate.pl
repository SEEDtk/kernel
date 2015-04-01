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
use ScriptUtils;
use Data::Dumper;
use SeedUtils;

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '',    ScriptUtils::ih_options(),
    [ 'reference|r=s','enome used for calibration', { required => 1 }]
);
my $r = $opt->reference;
print STDERR "using genome $r\n";
my $ih = ScriptUtils::IH($opt->input);
$/ = "\n//\n";
my @tuples;
while (defined($_ = <$ih>))
{
    if ($_ =~ /^(\d+\.\d+)\t(\S[^\n]*)\n---\n(.*)\n\/\//s)
    {
	my($g,$gs,$hits) = ($1,$2,$3);
	foreach my $x (split(/\n/,$hits))
	{
	    if ($x =~ /^(\d+)\t(\d+\.\d+)\t(\S.*\S)/)
	    {
		my($n,$poss_r,$rs) = ($1,$2,$3);
		if ($poss_r eq $r)
		{
		    push(@tuples,[$g,$n]);
		}
	    }
	}
    }
}
@tuples = sort { $b->[1] <=> $a->[1] } @tuples;
my @todo;
my $hi = 100000;
while (@tuples > 0)
{
    my $tuple;
    while (($tuple = shift @tuples) && ($hi < $tuple->[1])) {}
    if ($tuple)
    {
	my($g,$n)  = @$tuple;
	print STDERR "processing $n\t$g\n";
	my $dataD  = "/Users/rossoverbeek/Proj/SEEDtk/Data/Calibrate";
	$hi        = $n - 500;
	if ($n > 1000)
	{
	    &SeedUtils::run("perl fetch_seed_gto.pl -s pseed -g $g > $dataD/$g.gto");
	    &SeedUtils::run("perl fast_project_ref_to_GTO.pl -r $r -i $dataD/$g.gto > $dataD/out.$r.$g 2> $dataD/$r.$g.$n.stderr");
	}
    }
}

   
