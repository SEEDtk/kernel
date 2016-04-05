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
use SeedUtils;
use Data::Dumper;

=head1 Test Projection of a Set of Subsystems

    run_set_of_solid_subsystems -d DataDir < file-of-subsystem-names

This script evaluates the quality of projection for a set of subsystems on the coreSEED.

=head2 Output

The script builds subdirectories of the DataDir, each named by subsystem-id


=cut

my $opt = ScriptUtils::Opts('',
          ScriptUtils::ih_options(),
	  [ 'ksize|k=i','Kmer Size', { default => 8 }],
          [ 'data|d=s', 'Name or Data Directory to contain results', { required => 1 } ]
    );
my $ih         = ScriptUtils::IH($opt->input);
my $dataD      = $opt->data();
my $k          = $opt->ksize();
if (! -d $dataD) { mkdir($dataD,0777) }
if (! -d $dataD) { die "coud not make the Data Directory $dataD" }

my $ss;
while (defined($ss = <$ih>))
{
    chomp $ss;
    if (-d "$dataD/$ss")
    {
	print STDERR "$dataD/$ss already exists; skipping it\n";
    }
    else
    {
	mkdir("$dataD/$ss",0777) || die "could not make $dataD/$ss";
	&SeedUtils::run("perl eval_subsys_proj.pl -s $ss -e $dataD/$ss/errors -A $dataD/$ss/all.fasta --cv 1 -a kmer_annotation.pl 2> $dataD/$ss/stderr > $dataD/$ss/stdout");
    }
}
