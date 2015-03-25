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
#   fetch_seed_gto -s seed-to-get-it-from -g genome_ID > gto.file
#
use strict;
use warnings;
use JSON::XS;
use FIG_Config;
use ScriptUtils;
use Data::Dumper;
use GenomeTypeObject;

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
    '', 
    ScriptUtils::ih_options(),
    [ 'genome|g=s','ID of the genome', { required => 1 }],
    [ 'seed|s=s','name of the seed', { required => 1 }]
);
my $genome       = $opt->genome;
my $seed         = $opt->seed;
my $ih = ScriptUtils::IH($opt->input);

use LWP::Simple;
my $ref_text;

if ($seed eq "core")
{
    $ref_text = &LWP::Simple::get("http://core.theseed.org/FIG/genome_object.cgi?genome=$genome");
}
elsif ($seed eq "pseed")
{
    $ref_text = &LWP::Simple::get("http://pseed.theseed.org/genome_object.cgi?genome=$genome");
}

if ($ref_text)
{
    print $ref_text;
}
