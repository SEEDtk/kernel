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
use Shrub::Subsystems;
use SeedUtils;
use ServicesUtils;

=head1 Project Subsystems onto Genomes

    project_subsystems.pl [ options ] 

Project subsystems onto one or more genomes. The genomes can be input as a set of genome IDs or as a single
L<GenomeTypeObject> in JSON format.

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item json

If specified, the input is presumed to be a L<GenomeTypeObject> in JSON format rather than a tab-delimited file
with genome IDs.

=item col

The index (1-based) of the input column containing genome IDs. A value of C<0> indicates the last column. The
default is C<0>.

=item priv

Privilege level for functional assignments. The default is C<0>.

=back

=head2 Output Format

The output will consist of the input line (genome ID if the C<--json> option is specified) followed by a subsystem ID,
a variant code, a role, and a feature ID. Since there will be many subsystems per genome and many features per subsystem,
there will be many more output rows than input rows.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(), ScriptUtils::ih_options(),
        ['json', 'input is JSON-form GenomeTypeObject'],
        ['col|c=i', 'index (1-based) of input column', { default => 0 }],
        ['priv|p=i', 'functional assignment privilege level', { default => 0 }]
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# This hash will hold the projection data for each genome.
my %results;
# This list will contain [genomeID, inputLine] couplets.
my @inputs; 
# Determine our course of action.
if ($opt->json) {
    # Here we have a GenomeTypeObject. We read using SeedUtils in case it's a kbase one, which is
    # incompatible. The projection code is smart enough to handle both, but the GenomeTypeObject code
    # can't be so flexible.
    my $gto = SeedUtils::read_encoded_object($ih);
    # Simulate getting the genome ID from the input file.
    my $genomeID = ServicesUtils::json_field($gto, 'id');
    push @inputs, [$genomeID, [$genomeID]];
    # Compute the projection.
    $results{$genomeID} = Shrub::Subsystems::ProjectForGto($shrub, $gto);
} else {
    # Here we have a tab-delimited input file. Get the privilege level.
    my $priv = $opt->priv;
    # Loop through the input.
    while (! eof $ih) {
        # Get a batch of genomes.
        my @couplets = ScriptUtils::get_couplets($ih, $opt->col, 10);
        push @inputs, @couplets;
        # Compute their subsystems.
        for my $couplet (@couplets) {
            my $genomeID = $couplet->[0];
            $results{$genomeID} = Shrub::Subsystems::ProjectForGenome($shrub, $genomeID, $priv);
        }
    }
}
# Output the results.
for my $input (@inputs) {
    my ($genomeID, $line) = @$input;
    my $subData = $results{$genomeID};
    if ($subData) {
        for my $sub (sort keys %$subData) {
            my $projectionData = $subData->{$sub};
            my ($variant, $subRow) = @$projectionData;
            for my $subCell (@$subRow) {
                my ($role, $fid) = @$subCell;
                print join("\t", @$line, $sub, $variant, $role, $fid) . "\n";
            }
        }
    }
}