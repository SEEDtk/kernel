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

=head1 Find Bad Subsystem Variants

    find_bad_variants [ options ]

This script creates the bad-variants list for L<solid_projection_pipeline.pl>. The main class of bad variant is
one with suspiciously few roles. To find such variants, the script loops through the subsystems, looking for
subsystem rows than contain fewer cells than all other rows. If all of the rows for a variant are among these
small-cell rows, then we predict that something is wrong with the variant. An additional option includes variants
with suspicious names, such as C<dirty> or C<missing>, or names that end in C<.x>.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item maxbadrows

The maximum number of rows that can be considered interesting.  When we build our set, we look for rows that
contain fewer cells that all other rows. The value of this parameter is the maximum number of rows we will
look at. Essentially, we sort the rows by the number of filled cells and look at the set determined by the
number of rows indicated in this parameter. If no rows outside the set have the same number of cells, then
we consider the row problematic. A value of C<1> here (the default) means we only consider a row problematic
if it has fewer filled cells than all other rows.

=item variantsOnly

If specified, the bad variants will be listed, one per line. Otherwise, the output will include a list of the
identified genomes.

=item suspicious

If specified, variants with suspicious variant codes will also be listed.

=back

=head2 Output

The output is to the standard output, and is a single tab-delimited file. There are two types of records.
A I<variant record> represents a suspicious variant, and consists of three columns, (0) the subsystem ID,
(1) the subsystem name, and (2) the variant code. Following each variant record is one or more I<genome records>
listing the genomes that have subsystem rows in the variant. Each variant record begins with an empty column,
followed by (1) the genome ID and (2) the genome name.

=cut


my $startTime = time;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
                Shrub::script_options(),
                ['maxbadrows|m=i', 'maximum number of subsystem rows to consider as problematic', { default => 1 }],
                ['variantsOnly', 'only list the variants, not the genomes within'],
                ['suspicious', 'list variants with suspicious variant codes']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
my $stats = Stats->new();
# Get the command-line parameters.
my $max_bad_rows = $opt->maxbadrows;
my $variantsOnly = $opt->variantsonly;
my $suspicious = $opt->suspicious;
# Loop through the subsystems.
my %subs = map { $_->[0] => $_->[1] } $shrub->GetAll('Subsystem', '', [], 'id name');
my $totSS = scalar keys %subs;
for my $ss (sort keys %subs) {
    my $count = $stats->Add(subsystem => 1);
    print STDERR "Processing $ss ($count of $totSS).\n";
    # Get all the row information for this subsystem.
    my @rowData = $shrub->GetAll("Subsystem SubsystemRow Genome AND SubsystemRow Row2Cell Cell2Feature",
                          'Subsystem(id) = ? AND SubsystemRow(needs-curation) = ? ORDER BY SubsystemRow(variant-code), Cell2Feature(from-link)',
                          [$ss, 0], [qw(SubsystemRow(variant-code) SubsystemRow(id) Cell2Feature(from-link)
                          Cell2Feature(to-link) Genome(id) Genome(name))]);
    # This hash will count the number of filled cells for each subsystem row.
    my %rowSize;
    # This hash will track the variant code for each row.
    my %rowVC;
    # This hash will list the rows for each variant code.
    my %vcRow;
    # This hash will list the variant codes output as bad.
    my %vcOut;
    # This hash will contain the genome names and IDs for each row.
    my %rowGenome;
    # Loop through the row data.
    my $oldCell = '';
    for my $rowDatum (@rowData) {
        my ($vc, $row, $cell, $fid, $g, $gName) = @$rowDatum;
        $stats->Add(rowData => 1);
        # Only process if this is a new cell.
        if ($cell eq $oldCell) {
            $stats->Add(doubledCells => 1);
        } else {
            # Connect the variant code to the row.
            $rowVC{$row} = $vc;
            $vcRow{$vc}{$row} = 1;
            $rowGenome{$row} = [$g, $gName];
            # Count this cell.
            $rowSize{$row}++;
            $stats->Add(countedCells => 1);
            # Set up for the next iteration.
            $oldCell = $cell;
        }
    }
    # Sort the cells.
    my @rows = sort { $rowSize{$a} <=> $rowSize{$b} } keys %rowSize;
    # Check to see if we have enough rows to make the check plausible.
    if (scalar @rows <= $max_bad_rows) {
        print STDERR "Subsystem $ss is too small for consideration.\n";
        $stats->Add(subsysTooSmall => 1);
    } else {
        # Get the minimum number of cells per row for a row to be considered good.
        my $goodCount = $rowSize{$rows[$max_bad_rows]};
        # This hash will contain the variant codes that could be problematic.
        my %vcBad;
        for (my $i = 0; $i < $max_bad_rows && $rowSize{$rows[$i]} < $goodCount; $i++) {
            # Here we have a problematic row. It has too few cells.
            $vcBad{$rowVC{$rows[$i]}} = 1;
            $stats->Add(problemRow => 1);
        }
        # Loop through the variant codes found and verify they are unique to the small rows.
        for my $vc (keys %vcBad) {
            # Look for a row with this variant code that is outside the small set.
            my $vcGood;
            for my $row (keys %{$vcRow{$vc}}) {
                if ($rowSize{$row} >= $goodCount) {
                    $vcGood = 1;
                }
            }
            # If this variant code is NOT good, write it out.
            if ($vcGood) {
                $stats->Add(problemVariantPassed => 1);
            } else {
                $vcOut{$vc} = 1;
                $stats->Add(problemVariantFound => 1);
                PrintVariant($ss, $vc, \%rowGenome, $vcRow{$vc});
            }
        }
    }
    # Here we process suspicious variants that weren't found by the above
    # search.
    if ($suspicious) {
        for my $vc (keys %vcRow) {
            # Only proceed if we haven't put this code out before.
            if (! $vcOut{$vc}) {
                $stats->Add(variantChecked => 1);
                # Check for funny variant codes.
                if ($vc =~ /\.x$/i || $vc eq 'dirty' || $vc eq 'missing') {
                    $stats->Add(variantSuspicious => 1);
                    PrintVariant($ss, $vc, \%rowGenome, $vcRow{$vc});
                }
            }
        }
    }
}
$stats->Add(totalTime => (time - $startTime));
print STDERR "All done:\n" . $stats->Show();


sub PrintVariant {
    my ($ss, $vc, $rowGenome, $vcRow) = @_;
    print "$ss\t$subs{$ss}\t$vc\n";
    if (! $variantsOnly) {
        # Assemble the genome data for this variant code.
        my @genomes = map { $rowGenome->{$_} } keys %$vcRow;
        for my $genome (@genomes) {
            print join("\t", "", @$genome) . "\n";
        }
    }
}