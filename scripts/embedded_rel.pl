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
use ScriptUtils;
use Data::Dumper;
use File::Slurp;
use ScriptUtils;
use ScriptThing;

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '',
                     [], ScriptUtils::ih_options(),
		     ['embedded_file|f=s', { required => 1 }],
		     [ 'column|c=i', {} ]          ### How should I set default = last? (will -1 work/)
    );
my $ih = ScriptUtils::IH( $opt->input );
my $column = $opt->column;
my $embedded_file = $opt->embedded_file;

open(E,">$embedded_file") || die "could not open $embedded_file";
my $nxt_emb = 1;
while (defined($_ = <STDIN>))
{
    chomp;
    my @flds = split(/\t/,$_);
    my $rel = $flds[$column-1];
    $flds[$column-1] = $nxt_emb;
    print join("\t",@flds),"\n";
    my @e = split(/,/,$rel);
    my $y = 1;
    foreach my $x (@e)
    {
	print E join("\t",($nxt_emb,$y,$x)),"\n";
	$y++;
    }
    $nxt_emb++;
}
close(E);
