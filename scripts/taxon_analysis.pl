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

=head1 Analysis Role Frequency by Taxonomy

    taxon_analysis.pl [ options ]

This script determines the universal roles for taxonomic groupings. A I<universal role> is one that occurs singly 95% of the time in genomes
of a given taxonomic grouping. The script takes as input a matrix of the singly-occurring roles in each known good genome. This matrix can
be used to compute the universal roles for each taxonomic grouping in the Shrub.

Progress messages are sent to STDERR. The output will be a tab-delimited file, one row per taxonomic grouping, containing the taxonomy ID,
the group name, and the universal roles found.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> (to specify the database) and L<ScriptUtils/ih_options> (to specify the
input file) plus the following.

=over 4

=item min

Minimum percentage of genomes in a taxonomic grouping that must contain a role for it to be considered universal. The default is C<95>.

=item size

Minimum number of genomes in a taxonomic grouping required for the grouping to be considered worthwhile. The default is C<100>.

=item rMin

Minimum number of roles for a taxonomic grouping for it to be considered useful in determining completeness.

=back

The input file contains a header row with the role IDs in it. The first column of each data line should contain the genome ID, the second the
taxonomic ID, and the remaining columns a C<1> for a singly-occurring role and C<0> otherwise.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['min|m=f', 'minimum threshold percentage', { default => 95 }],
        ['size|s=i', 'minimum number of genomes per group', { default => 100 }],
        ['rMin|r=i', 'minimum number of roles for a useful group', { default => 100 }]
        );
my $stats = Stats->new();
# Get the merge file. We need to map each old taxon to its new version.
open(my $mh, "<$FIG_Config::data/Inputs/Others/merged.dmp") || die "Could not open merged.dmp: $!";
my %merge;
while (! eof $mh) {
    my $line = <$mh>;
    if ($line =~ /(\d+)\t\|\t(\d+)/) {
        $merge{$1} = $2;
        $stats->Add(taxonMergeRead => 1);
    }
}
# Connect to the database.
print STDERR "Connecting to the database.\n";
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Read the header line.
my (undef, undef, @roles) = ScriptUtils::get_line($ih);
# This hash counts the genomes in each grouping.
my %taxCounts;
# This hash contains the taxonomic grouping names.
my %taxNames;
# This hash contains the array of role counts for each group.
my %taxRoles;
# This hash maps each taxonomic grouping to its parent.
my %taxParent;
# Use this to time trace messages.
my $count = 0;
print STDERR "Processing genomes from input.\n";
# Loop through the input.
while (! eof $ih) {
    my ($genome, $taxon, @array) = ScriptUtils::get_line($ih);
    while ($taxon) {
        # Make sure we have the latest version of this tax ID.
        if ($merge{$taxon}) {
            $taxon = $merge{$taxon};
            $stats->Add(taxonMapped => 1);
        }
        # Get the taxon name.
        my $taxName = $taxNames{$taxon};
        if (! $taxName) {
            # This is a new group. Get its parent and name.
            my ($taxInfo) = $shrub->GetAll('TaxonomicGrouping IsInTaxonomicGroup', 'TaxonomicGrouping(id) = ?', [$taxon], 'scientific-name domain IsInTaxonomicGroup(to-link)');
            if (! $taxInfo) {
                print STDERR "Taxonomic grouping $taxon not found for $genome.\n";
            } else {
                my ($name, $domain, $parent) = @$taxInfo;
                $taxNames{$taxon} = $name;
                $taxParent{$taxon} = ($domain ? '' : $parent);
                $taxName = $name;
                $stats->Add(taxGroupFound => 1);
            }
        }
        # Only proceed if this is a valid grouping.
        if (! $taxName) {
            $stats->Add(taxonInvalid => 1);
            $taxon = '';
        } else {
            # Add us into this group.
            if ($taxRoles{$taxon}) {
                for (my $i = 0; $i < scalar @roles; $i++) { $taxRoles{$taxon}[$i] += $array[$i] }
                $stats->Add(taxGroupReused => 1);
            } else {
                $taxRoles{$taxon} = [ @array ];
                $stats->Add(taxGroup => 1);
            }
            $taxCounts{$taxon}++;
            # Get the next group.
            $taxon = $taxParent{$taxon};
        }
    }
    $stats->Add(genomeProcessed => 1);
    $count++;
    print STDERR "$count genomes processed.\n" if ($count % 1000 == 0);
}
# Now we run through the taxonomic groupings producing the output.
print STDERR scalar(keys %taxRoles) . " total groups found.\n";
my @goodTaxes = sort grep { $taxCounts{$_} >= $opt->size } keys %taxRoles;
print STDERR "Producing output.  " . scalar(@goodTaxes) . " sufficiently large groups found.\n";
print join("\t", 'TaxID', 'Name', @roles) . "\n";
for my $taxon (@goodTaxes) {
    $stats->Add(checkedGroups => 1);
    my $gTotal = $taxCounts{$taxon};
    my $gMin = int(($opt->min * $gTotal + 99) / 100);
    my @output = ($taxon, $taxNames{$taxon});
    # This counts the marker roles.
    my $roleCount = 0;
    # Loop through all the roles, putting 0 normally and 1 for markers.
    for my $val (@{$taxRoles{$taxon}}) {
        my $mark = (($val >= $gMin) ? 1 : 0);
        $roleCount += $mark;
        push @output, $mark;
    }
    if ($roleCount >= $opt->rmin) {
        print join("\t", @output) . "\n";
        $stats->Add(goodGroups => 1);
        $stats->Add(markers => $roleCount);
    } else {
        $stats->Add(badGroups => 1);
    }
}
print STDERR "All done.\n" . $stats->Show();
