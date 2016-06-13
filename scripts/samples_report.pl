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

=head1 Output Sample Bin Analysis

    samples_report.pl [ options ] site1 site2 ...

This script lists all the metagenomic samples in the database for a given list of sites, one line per bin.
For the bin, it will display the genus of the bin and the percentage of the total content the bin comprises.

The output is a tab-delimited file, each record containing (0) a sample ID, (1) the name of the current bin's genus,
(2) the percent of the sample's binned contigs in the current bin, and (3) the site ID.

=head2 Parameters

The positional parameters are the names of the sites whose bins are to be listed.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item project

If specified, the name of a project. Only sites from the specified project will be included.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        ['project=s', 'project to which the queries should be restricted'],
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the project ID.
my $project = $opt->project;
# Get the site list.
my @sites = @ARGV;
# Loop through the sites.
for my $site (@sites) {
    my @bins = $shrub->GetAll('Site2Sample Sample2Bin Bin Bin2Taxonomy TaxonomicGrouping',
            'Site2Sample(from-link) = ?', [$site], 'Sample2Bin(from-link) Bin(id) Bin(dna-size) TaxonomicGrouping(scientific-name)');
    # Get the DNA length for each sample and throw out bins from the wrong project.
    my %binLength;
    my @actualBins;
    my $projLen;
    if ($project) {
        $project .= '.';
        $projLen = length $project;
    } 
    for my $bin (@bins) {
        my ($sample, undef, $len) = @$bin;
        if (! $project || substr($sample, 0, $projLen) eq $project) {
            $binLength{$sample} += $len;
            push @actualBins, $bin;
        } 
    }
    # Now output the bins. We sort by sample ID followed by dna size descending
    for my $bin (sort { ($a->[0] cmp $b->[0]) || ($b->[2] <=> $a->[2]) } @bins) {
        my ($sample, $binID, $dnaSize, $genus) = @$bin;
        my $percent = $dnaSize * 100 / $binLength{$sample};
        print "$sample\t$genus\t$percent\t$site\n";
    }
}