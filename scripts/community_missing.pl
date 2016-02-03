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
use SampleContig;


=head1 Describe Script Here

    community_missing.pl [ options ] binFile nodeFile

This script looks at the bin output file from a L<community_pipeline.pl> run and determines which nodes are missing
and what their lengths and coverages are. The coverages are rounded to the nearest whole number and the lengths
rounded to the nearest 1000.

=head2 Parameters

There are two positional parameters-- the bin output file and the original FASTA file containing the nodes.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binFile nodeFile');
# Create the statistics object.
my $stats = Stats->new();
# Get the parameters.
my ($binFile, $nodeFile) = @ARGV;
if (! $binFile) {
    die "No bins input file specified.";
} elsif (! $nodeFile) {
    die "No nodes FASTA file specified.";
}
# Read the bin file and get a hash of the nodes placed.
my %binned;
print STDERR "Reading bin file $binFile.\n";
open(my $bh, "<", $binFile) || die "Could not open bin file: $!";
while (! eof $bh) {
    # Get a node ID.
    my $line = <$bh>;
    chomp $line;
    if ($line ne "//") {
        $binned{$line} = 1;
        # Skip over the node information.
        while (! eof $bh && $line) {
            $line = <$bh>;
            chomp $line;
        }
    }
}
close $bh;
# Now read the node file and analyze the nodes not found.
print STDERR "Reading node file $nodeFile.\n";
open(my $nh, "<", $nodeFile) || die "Could not open node file: $!";
while (! eof $nh) {
    my $line = <$nh>;
    if ($line =~ /^>(\S+)/) {
        my $nodeID = $1;
        my $contig = SampleContig->new($nodeID);
        my $covg = int(($contig->covg() + 5) / 10) . "0";
        my $size = int(($contig->len() + 5000) / 10000) . "0K";
        if ($binned{$nodeID}) {
            $stats->Add(foundBin => 1);
            $stats->Add("foundCovg$covg" => 1);
            $stats->Add("foundSize$size" => 1);
        } else {
            $stats->Add(skippedBin => 1);
            $stats->Add("skippedCovg$covg" => 1);
            $stats->Add("skippedSize$size" => 1);
        }
    }
}
print "All done.\n" . $stats->Show();
