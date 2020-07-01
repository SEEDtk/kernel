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
use Stats;
use Bin;

=head1 Analyze Bin Improvement Contigs

    analyze_bin_contigs.pl [ options ] binDir

This script will read the bins.json files from a bin directory and list the contigs associated with each bin.  For each contig, its ID and length will be output, along with
an indication of whether or not the contig was removed during improvement.

=head2 Parameters

The positional parameter is the name of the binning directory.

The command-line options are the following

=over 4

=item verbose

Display progress information on STDERR.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ["verbose|debug|v", "display progress on STDERR"]
        );
my $debug = $opt->verbose;
# Get the binning directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -s "$binDir/bins.rast.json") {
    die "$binDir does not have a bins.rast.json file.";
}
# Read in the final bins.
print STDERR "Reading bins.rast.json for $binDir.\n" if $debug;
my $binsRast = binHash(Bin::ReadBins("$binDir/bins.rast.json"));
# Read in the initial bins.
print STDERR "Reading bins.json for $binDir.\n" if $debug;
my $binsInit = binHash(Bin::ReadBins("$binDir/bins.json"));
# Write column headers.
print "sample\tbin_name\tcontig_id\tlen\tcovg\tfinal\n";
# The "binHash" method organized each bin list into a hash.  We loop through the hashes in parallel.
for my $binTaxId (sort keys %$binsRast) {
    my $initBin = $binsInit->{$binTaxId};
    my $rastBin = $binsRast->{$binTaxId};
    my $name = $rastBin->name;
    if (! $initBin) {
        print STDERR "WARNING: no initial version found for $binTaxId: $name.\n" if $debug;
    } else {
        print STDERR "Processing bin $name\n" if $debug;
        # Create a hash of the final-bin contigs.
        my %rastContigs = map { $_ => 1 } $rastBin->contigs;
        # Loop through the initial-bin contigs.
        my $initContigs = $initBin->all_contigs;
        my @sorted = sort { $b->len <=> $a->len } @$initContigs;
        for my $contig (@sorted) {
            my $contigId = $contig->id;
            my $contigLen = $contig->len;
            my $contigCovg = $contig->covg;
            my $kept = ($rastContigs{$contigId} ? "X" : "");
            print "$binDir\t$name\t$contigId\t$contigLen\t$contigCovg\t$kept\n";
        }
    }
}
print STDERR "All done.\n" if $debug;


## Hash the bins by taxon ID.
sub binHash {
    my ($binList) = @_;
    my $retVal = {};
    for my $bin (@$binList) {
        $retVal->{$bin->{taxonID}} = $bin;
    }
    return $retVal;
}
