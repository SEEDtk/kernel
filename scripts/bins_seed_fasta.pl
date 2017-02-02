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
use Shrub;
use ScriptUtils;

=head1 Generate Seed Protein FASTA for Binning

    bins_seed_fasta.pl [ options ] role

This script generates a FASTA file that can be used to prime the binning process. It outputs a FASTA
file containing protein sequences from specific genomes for the named universal protein.

=head2 Parameters

The positional parameter is the name of the protein to use.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item genomes

A comma-delimited list of the genomes from which to extract the protein sequences.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('role',
                Shrub::script_options(),
                ['seedgenome|G=s', 'ID of the genome to seed the bins', { default => '83333.1,36870.1,224308.1,1148.1,64091.1,69014.3,83332.1,115711.7,187420.1,224326.1,243273.1,4932.3' }],
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the protein ID and the genome list.
my ($prot) = @ARGV;
my @genomes = split /,/, $opt->seedgenome;
# Retrieve the proteins.
my $filter = "Function2Feature(security) = ? AND Function2Feature(from-link) = ? AND Feature2Genome(to-link) IN (" .
        join(', ', map { '?' } @genomes) . ")";
my $parms = [0, $prot, @genomes];
my @protPairs = $shrub->GetAll('Function2Feature Feature Feature2Genome AND Feature Protein', $filter, $parms,
        'Feature(id) Protein(sequence)');
for my $protPair (@protPairs) {
    my ($id, $sequence) = @$protPair;
    print ">$id\n$sequence\n";
}
