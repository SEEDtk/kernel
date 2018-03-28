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

=head1 Analyze Role Frequency by Taxonomy

    taxon_analysis.pl [ options ] outDir

This script determines the universal roles for taxonomic groupings. A I<universal role> is one that occurs singly 95% of the time in genomes
of a given taxonomic grouping. The script takes as input a matrix of the singly-occurring roles in each known good genome (the output of
L<p3-role-matrix.pl>). This matrix can be used to compute the universal roles for each taxonomic grouping in the Shrub.

Progress messages are sent to the standard output. The output directory will contain the following files, all tab-delimited.

=over 4

=item taxon.tbl

A list of the useful taxonomic groupings, one per line, containing the taxonomic ID, the name, and an indicator (C<1> or C<0>) of the
required roles, one role per column.

=item sizes.tbl

A list of all taxonomic groupings, one per line, containing (0) the taxonomic ID, (1) the number of genomes, and (2) the number of
required roles.

=item roles.tbl

A list of the useful taxonomic groupings, one per line, containing the taxonomic ID, the name, and a list of required role IDs.

=back

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

=item merge

The name of the NCBI taxonomy database file containing mappings from old taxon IDs to new ones.

=back

The input file contains a header row with the role IDs in it. The first column of each data line should contain the genome ID, the second the
taxonomic ID, and the remaining columns a C<1> for a singly-occurring role and C<0> otherwise.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['min|m=f', 'minimum threshold percentage', { default => 95 }],
        ['size|s=i', 'minimum number of genomes per group', { default => 100 }],
        ['rMin|r=i', 'minimum number of roles for a useful group', { default => 100 }],
        ['merge=s', 'NCBI taxonomy merge file', { default => "$FIG_Config::data/Inputs/Other/merged.dmp" }]
        );
my $stats = Stats->new();
# Verify the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir);
}
# Get the merge file. We need to map each old taxon to its new version.
open(my $mh, '<', $opt->merge) || die "Could not open merged.dmp: $!";
my %merge;
while (! eof $mh) {
    my $line = <$mh>;
    if ($line =~ /(\d+)\t\|\t(\d+)/) {
        $merge{$1} = $2;
        $stats->Add(taxonMergeRead => 1);
    }
}
# Connect to the database.
print "Connecting to the database.\n";
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
print "Processing genomes from input.\n";
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
    print "$count genomes processed.\n" if ($count % 100 == 0);
}
# Now we run through the taxonomic groupings producing the output. Create the output files.
print scalar(keys %taxRoles) . " total groups found.\n";
open(my $sh, ">$outDir/sizes.tbl") || die "Could not open sizes.tbl: $!";
print $sh join("\t", 'TaxID', 'Name', 'Size', '#Roles') . "\n";
open(my $rh, ">$outDir/roles.tbl") || die "Could not open roles.tbl: $!";
print $rh join("\t", 'TaxID', 'Name', 'Size', 'Required Roles') . "\n";
open(my $th, ">$outDir/taxon.tbl") || die "Could not open taxon.tbl: $!";
print $th join("\t", 'TaxID', 'Name', @roles) . "\n";
# Get the number of roles.
my $nRoles = scalar @roles;
for my $taxon (sort keys %taxRoles) {
    $stats->Add(checkedGroups => 1);
    my $gTotal = $taxCounts{$taxon};
    my $gMin = int(($opt->min * $gTotal + 99) / 100);
    my $tCounts = $taxRoles{$taxon};
    my @header = ($taxon, $taxNames{$taxon});
    # This will be the <1,0> list.
    my @marks;
    # This will be the role ID list.
    my @ids;
    # This counts the marker roles.
    my $roleCount = 0;
    # Loop through all the roles for this taxon.
    for (my $i = 0; $i < $nRoles; $i++) {
        my $val = $tCounts->[$i];
        my $role = $roles[$i];
        my $mark = (($val >= $gMin) ? 1 : 0);
        push @marks, $mark;
        if ($mark) {
            $roleCount++;
            push @ids, $role;
        }
    }
    print $sh join("\t", @header, $taxCounts{$taxon}, $roleCount) . "\n";
    if ($taxCounts{$taxon} < $opt->size) {
        $stats->Add(groupTooSmall => 1);
        print "$taxon: $taxNames{$taxon} is too small ($taxCounts{$taxon} genomes).\n";
    } elsif ($roleCount < $opt->rmin) {
        $stats->Add(groupBad => 1);
        print "$taxon: $taxNames{$taxon} has too few marker roles ($roleCount).\n";
    } else {
        print "$taxon: $taxNames{$taxon} is good-- $roleCount markers found in $taxCounts{$taxon} genomes.\n";
        print $th join("\t", @header, @marks) . "\n";
        print $rh join("\t", @header, $taxCounts{$taxon}, @ids) . "\n";
        $stats->Add(groupsGood => 1);
        $stats->Add(markers => $roleCount);
    }
}
print "All done.\n" . $stats->Show();
