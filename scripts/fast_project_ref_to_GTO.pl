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

    fast_project_ref_to_GTO -r RererenceGenomeId [-k kmer-sz] < gto [ options ]

project a reference genome to call features in a new genome (the gto)


=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -r ReferenceGenomeId

a path to a SEED genome for the reference genome

=back

the gto must include the genome ID and the genetic code

=cut


# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
	'',
          ScriptUtils::ih_options(),
	[ 'reference|r=s', 'Id of the Reference Genome' ],
	[ 'kmersize|k=i',  'Number of base pairs per kmer', { default => 30 }]
);
my $ih       = ScriptUtils::IH($opt->input);
my $g_gto    = &GenomeTypeObject::create_from_file($ih);

my $refId    = $opt->reference;
use LWP::Simple;

my $ref_text = &LWP::Simple::get("http://core.theseed.org/FIG/genome_object.cgi?genome=$refId");
my $ref_gto  = JSON::XS->new;
$ref_gto->decode($ref_text);
bless($ref_gto,'GenomeTypeObject');

my $k        = $opt->kmersize;

my $genetic_code = $g_gto->{genetic_code};

my @ref_tuples = map { [$_->{contig_id},'',$_->{dna}] } @{$ref_gto->contigs};
my @g_tuples   = map { [$_->{contig_id},'',$_->{dna}] } @{$g_gto->contigs};

my $map = &MapToRef::build_mapping($k, \@ref_tuples, \@g_tuples );
&MapToRef::build_features( $map, my $refD, my $genomeD, \@g_tuples, $genetic_code );

