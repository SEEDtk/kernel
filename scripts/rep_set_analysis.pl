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
use Stats;

=head1 Analyze a Set of Representative Genomes

    rep_set_analysis.pl [ options ]

This script analyzes the output from a representative genome run and sorts them into Core, non-Core Shrub, and
PATRIC genomes. This information is then used to modify the Shrub load.

The input file has two header lines followed by tab-delimited genome lines. The first column is the genome index,
the second is the genome ID, and the third is the genome name.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> and L<ScriptUtils/ih_options> (to indicate the rep_server
output).

=cut

my $stats = Stats->new();
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        ScriptUtils::ih_options());
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# These hashes store the different types of genomes.
my %genomes = (core => {}, shrub => {}, missing => {});
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(lineIn => 1);
    if ($line =~ /^\d+\t(\d+\.\d+)\t(.+)/) {
        my ($id, $name) = ($1, $2);
        $stats->Add(genomes => 1);
        # Determine the genome type.
        my ($core) = $shrub->GetFlat('Genome', 'Genome(id) = ?', [$id], 'core');
        my $type;
        if (! defined $core) {
            # Not in Shrub.
            $type = 'missing';
        } elsif ($core) {
            $type = 'core';
        } else {
            $type = 'shrub';
        }
        # Count this genome.
        $stats->Add($type => 1);
        $genomes{$type}{$id} = $name;
    }
}
# Now produce the output, listing the genomes of each type.
for my $type (sort keys %genomes) {
    my $gHash = $genomes{$type};
    print uc($type) . " Genomes\n";
    for my $genome (sort keys %$gHash) {
        print "$genome\t$gHash->{$genome}\n";
    }
    print "\n";
}
print "Statistics:\n" . $stats->Show();
