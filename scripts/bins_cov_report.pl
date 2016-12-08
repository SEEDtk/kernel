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
use ScriptUtils;
use GenomeTypeObject;

=head1 Produce Binning Coverage Report

    bin_cov_report.pl [ options ] dir bin

This script processes a single sample and produces a report on the bins. For each bin it will list the name, the mean coverage
overall, and the mean coverage for the contigs containing seed proteins.

=head2 Parameters

There are two positional parameters-- the name of the directory containing the samples, and the ID of the sample to
be processed. If the sample ID is C<all>, then all samples will be processed.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir sample',
        );
# Get the bin directory.
my ($dir, $sampleID) = @ARGV;
if (! $sampleID) {
    die "No sample ID specified.";
} elsif (! $dir) {
    die "No samples directory specified.";
} elsif (! -d $dir) {
    die "Invalid samples directory $dir: $!";
} else {
    my @samples;
    if ($sampleID eq 'all') {
        opendir(my $dh, $dir) || die "Could not open sample directory: $!";
        @samples = sort grep { substr($_,0,1) ne '.' && -d "$dir/$_" } readdir $dh;
    } else {
        push @samples, $sampleID;
    }
    for my $sample (@samples) {
        my $sampleDir = "$dir/$sample";
        # Only process completed samples.
        if (-f "$sampleDir/bins.rast.json") {
            # Read in the contig IDs that represent seed protein hits.
            my %seedContigs;
            open(my $ih, '<', "$sampleDir/bins.found.tbl") || die "Could not open bins.found.tbl: $!";
            while (! eof $ih) {
                my $line = <$ih>;
                if ($line =~ /^(NODE_\S+)\t/) {
                    $seedContigs{$1} = 1;
                }
            }
            # Loop through the bins.
            my $i = 1;
            while (-f "$sampleDir/bin$i.gto") {
                # Read in this bin.
                my $gto = GenomeTypeObject->create_from_file("$sampleDir/bin$i.gto");
                # Get the name.
                my $name = $gto->{scientific_name};
                # Get the list of contigs.
                my $contigList = $gto->{contigs};
                # Loop through the contigs, accumulating averages.
                my ($tot, $count, $seedTot, $seedCount) = (0, 0, 0, 0);
                for my $contig (@$contigList) {
                    # Get this contig's ID and coverage.
                    my $contigID = $contig->{id};
                    if ($contigID =~ /cov_(.+)/) {
                        my $cov = $1;
                        $count++;
                        $tot += $cov;
                        if ($seedContigs{$contigID}) {
                            $seedTot += $cov;
                            $seedCount++;
                        }
                    }
                }
                # Compute the averages.
                my $avg = ($count ? $tot / $count : 0);
                my $seedAvg = ($seedCount ? $seedTot / $seedCount : 0);
                # Output this bin.
                print join("\t", $sample, $name, $avg, $seedAvg) . "\n";
                # Move to the next bin.
                $i++;
            }
        }
    }
}
