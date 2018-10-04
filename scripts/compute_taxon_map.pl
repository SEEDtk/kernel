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

=head1 Compute Taxon Map for Taxon/Role Table

    compute_taxon_map.pl [ options ]

This script takes an input file of taxonomic IDs. It outputs a file that maps all taxonomic IDs to the smallest containing taxonomic ID found in the
input file. The output file can then be used to map any represented taxonomic ID to the best match in the input file. Unrepresented taxonomies will not be included.
The basic strategy is to work our way up from each ID in the taxonomy tree and return the first group we find from the input.

=head2 Parameters

There are no positional parameters.

The input file will have headers and should have a taxonomic ID in its first column.

The command-line options are those found in L<Shrub/script_options> and L<ScriptUtils/ih_options> (to specify the input).

=head2 Output

The standard output will be a two-column file, the first column being an input taxonomic ID and the second being the appropriate representative from the input
file.

Progress information will be written to the standard error output.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        );
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# We want to create a hash of the input.
my %taxonsIn;
print STDERR "Reading input file.\n";
# Discard the header line.
my $line = <$ih>;
# Loop through the data lines.
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(lineIn => 1);
    if ($line =~ /^(\d+)\t([^\t]+)/) {
        $taxonsIn{$1} = $2;
        $stats->Add(taxonIn => 1);
    } else {
        die "Invalid input line: $line";
    }
}
print STDERR scalar(keys %taxonsIn) . " taxonomic group IDs found.\n";
# Output the header for the output file.
print join("\t", qw(taxonID group name)) . "\n";
# Now we loop through all the taxonomic groupings.
my $q = $shrub->Get('TaxonomicGrouping', '', [], 'id scientific-name');
print STDERR "Processing taxonomic groups from database.\n";
my ($count, $mapped) = (0, 0);
while (my $record = $q->Fetch()) {
    $stats->Add(dbTaxon => 1);
    my ($id, $name) = $record->Values('id scientific-name');
    # This will be set to the group found.
    my $found;
    my $nextID = $id;
    # Crawl up the tree.
    while ($nextID > 1 && ! $found) {
        if ($taxonsIn{$nextID}) {
            $found = $nextID;
            $stats->Add(groupFound => 1);
            $mapped++;
        } else {
            ($nextID) = $shrub->GetFlat('IsInTaxonomicGroup', 'IsInTaxonomicGroup(from-link) = ?', [$nextID], 'to-link');
            $stats->Add(groupChecked => 1);
        }
    }
    # Output the result.
    if (! $found) {
        $stats->Add(groupNotMapped => 1);
    } else {
        print join("\t", $id, $found, $taxonsIn{$found}) . "\n";
        $stats->Add(groupMapped => 1);
    }
    $count++;
    if ($count % 1000 == 0) {
        print STDERR "$count groups processed, $mapped mapped.\n";
    }
}
print STDERR "All done.\n" . $stats->Show();
