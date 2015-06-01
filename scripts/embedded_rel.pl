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

=head1 Extract an Embedded, comma-separated Set

A comma-separated embedded set is a set of values separated by commas, occurring as
a single value in a tab-separated table.  Thus,

    PossRnaDegrClus:125	fig|83333.1.peg.2702,fig|83333.1.peg.2703
    RiboProtZincRequ:208	fig|83333.1.peg.292,fig|83333.1.peg.293

might be two rows in a tab-sepataed table in which the second column has
embedded sets as values.

This utility creates two files: the original file in which the embedded 
sets are replaced with an integer (a unique value for each row of the table).
The embedded sets are expanded in a separate file as

         [integer-id-of-embedded-set,value]

Thus, the two rows above would be converted to

MAIN-TABLE

    PossRnaDegrClus:125	1
    RiboProtZincRequ:208	2

EMBEDDED SETS
    1	fig|83333.1.peg.2702
    1   fig|83333.1.peg.2703
    2	fig|83333.1.peg.292
    2	fig|83333.1.peg.2933

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

-c N  gives the column containing the embedded set

-f File  gives the file to which the embedded sets are written.

=cut

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
    foreach my $x (@e)
    {
	print E join("\t",($nxt_emb,$x)),"\n";
    }
    $nxt_emb++;
}
close(E);
