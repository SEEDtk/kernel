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
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use gjo::seqlib;
use MapToRef;

=head1 project reference genome to a close strain

    fast_project -r RererenceGenomeDir -g SkeletalGenomeDir [ options ]

project a reference genome to call features


=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -r ReferenceGenomeDir

a path to a SEED genome directory for the reference genome

=item -g SkeletalGenomeDir

a path to a skeletal SEED genome directory that must include

    contigs
    GENETIC_CODE (if not 11)

=back

=cut


# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
	'',
	[ 'reference|r=s', 'Path to Reference Genome Directory' ],
	[ 'genome|g=s',    'Path to a skeletal genome directory' ],
	[ 'kmersize|k=i',  'Number of base pairs per kmer', { default => 30 }]
);
my $refD    = $opt->reference;
my $genomeD = $opt->genome;
my $k       = $opt->kmersize;

my $genetic_code = &MapToRef::get_genetic_code($genomeD);

my @ref_tuples = &gjo::seqlib::read_fasta("$refD/contigs");
my @g_tuples   = &gjo::seqlib::read_fasta("$genomeD/contigs");

my $map = &MapToRef::build_mapping($k, \@ref_tuples, \@g_tuples );
&MapToRef::build_features( $map, $refD, $genomeD, \@g_tuples, $genetic_code );

