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

=head1 Compare Taxonomic Bin Expectations

    bins_custom_check.pl [ options ] expectFile

This script compares a bin report to a list of expected bins. Each expected bin is identified by a taxonomic grouping name
and ID. The actual bin matches if it is on the same chain as the expected bin. When comparing an actual and expected taxonomic
grouping, the I<distance> is the number of taxonomic levels between the expected group and the actual group on the chain. We want
the match with the smallest distance.

At the end, we list the matches, the unmatched expected bins, and the unmatched actual bins.

=head2 Parameters

The positional parameter is the name of the file containing the expectations.

The standard input file is a C<bins.report.txt> from a binning run.

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('expectFile',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# This will list the expected bins.
my @expected;
# This will be TRUE for each claimed bin.
my @claimed;
# This will map each taxonomic grouping ID to a list of [bin, distance] pairs.
my %taxonomyMap;
# Get the expectation file.
my ($expectFile) = @ARGV;
if (! $expectFile) {
    die "No expectation file specified.";
} elsif (! -s $expectFile) {
    die "Expectation file $expectFile not found or empty.";
} elsif (! open(my $eh, '<', $expectFile)) {
    die "Could not open expectation file: $!";
} else {
    while (! eof $eh) {
        my ($taxName, $taxID) = (<$eh> =~ /^(\S+)\s+(\d+)/);
        if (! $taxID) {
            die "Invalid format on expectation file line $..";
        } else {
            # Record this bin as unclaimed.
            push @expected, $taxName;
            push @claimed, 0;
            # Remember the bin number.
            my $binN = $#claimed;
            # Get the taxonomic hierarchy. We only go down, since the expectation is genus or higher, and the
            # actual is genus or species.
            my @taxons = ([$taxID, 0]);
            while (my $newTaxInfo = pop @taxons) {
                my ($taxID, $dist) = @$newTaxInfo;
                push @{$taxonomyMap{$taxID}}, [$binN, $dist];
                $dist++;
                my @taxIDs = $shrub->GetFlat('IsTaxonomicGroupOf', 'IsTaxonomicGroupOf(from-link) = ?', [$taxID], 'to-link');
                push @taxons, map { [$_, $dist] } @taxIDs;
            }
        }
    }
}
# Sort the map entries.
for my $taxID (keys %taxonomyMap) {
    my $possibles = $taxonomyMap{$taxID};
    my @sorted = sort { $a->[1] <=> $b->[1] } @$possibles;
    $taxonomyMap{$taxID} = \@sorted;
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the input. We need two lines from the report: the header line and the taxon line.
my ($binID, $binNode, $binQual);
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^BIN\s+(\d+).+?(NODE[^,]+).+?quality\s+(\d\.\d\d\d)/) {
        ($binID, $binNode, $binQual) = ($1, $2, $3);
    } elsif ($line =~ /NCBI\s+taxon\s+(\d+):\s+(.+)/) {
        my ($taxon, $name) = ($1, $2);
        # Now we have all we need to process this bin. Get the best possible match.
        # It cannot have already been claimed.
        my $possibles = $taxonomyMap{$taxon} // [];
        my $match = shift @$possibles;
        while ($match && $claimed[$match->[0]]) {
            $match = shift @$possibles;
        }
        if (! $match) {
            $match = "NO MATCH";
        } else {
            my $binM = $match->[0];
            $match = "$binM. $expected[$binM]";
            $claimed[$binM] = 1;
        }
        # Output the result.
        print "BIN $binID ($binNode)\t$binQual\t$taxon\t$name\t$match\n";
    }
}
# Now print the unclaimed bins.
print "\nUNCLAIMED BINS\n";
for (my $i = 0; $i < @claimed; $i++) {
    if (! $claimed[$i]) {
        print "$i. $expected[$i]\n";
    }
}
