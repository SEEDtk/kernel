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
use Bin;
use Stats;
use Cwd;

=head1 Analyze Clustering Results Against Bins

    cluster_bin_analysis.pl [ options ] binfile workDir

This script compares the B<bins.json> file from the L<bins_create.pl> script to a file of clustering
information for the same contigs (that is, the ones in the bins and perhaps others that had been rejected).
The clustering information should be in the form of output from L<svc_cluster_pairs.pl>, that is, a count
followed by a list of IDs in tab-delimited format. What we want to know is how many bins are represented
in each cluster and how many contigs are participating in the clusters.

=head2 Parameters

The first positional parameter is the name of the B<bins.json> file from L<bins_create,pl>. This tells us 
which contigs are in each bin. The second parameter is the name of the working directory. If it is omitted,
the current directory is assumed. The standard input should be the clustering file from L<svc_cluster_pairs.pl>.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item details

If specified, should be a number of bins. The details of any cluster with the specified number of bins or more
will be output to the C<details.tbl> file in the working directory. The file will be tab-delimited with each
record containing (0) a cluster ID, (1) a bin ID, and (2) the ID of a contig in the cluster that is in the
given bin.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binfile workDir', ScriptUtils::ih_options(),
        ['details=i', 'if specified, number of bins required to output bin details'],
        );
# Create the statistics object.
my $stats = Stats->new();
# Get the bin file and the working directory.
my ($binFile, $workDir) = @ARGV;
# Verify the working directory.
if (! $workDir) {
    $workDir = cwd();
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.";
}
# Read in the bin file.
my %binMap;
if (! $binFile) {
    die "A bins.json file is required.";
} elsif (! -f $binFile) {
    die "Invalid bin file name $binFile.";
} else {
    my $binList = Bin::ReadBins($binFile);
    # Create a map of contigs to bins.
    for my $bin (@$binList) {
       $stats->Add(binsFound => 1);
       # Get the list of contigs in this bin.
       my @contigs = $bin->contigs;
       # Only proceed if the bin is nontrivial.
       if (scalar @contigs > 1) {
           $stats->Add(bigBins => 1);
           # Put the contigs in the hash.
           my $binLabel = $bin->contig1;
           for my $contig (@contigs) {
               $binMap{$contig} = $binLabel;
               $stats->Add(binnedContigs => 1);
           }
       }
    }
}
# Check for the details option.
my $details = $opt->details;
my $oh;
if (defined $details) {
    open($oh, ">$workDir/details.tbl") || die "Could not open details file: $!";
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# This will be used as the cluster ID.
my $clusterID = 1;
# Loop through it.
while (! eof $ih) {
    # Get the current line and chop off the CR/LF stuff. (chomp does not always work, tragically)
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    # Extract the columns.
    my ($clusterCount, @contigs) = split /\t/, $line;
    $stats->Add(clustersIn => 1);
    # This will track the number of contigs in each bin.
    my %bins;
    # This will track the actual contigs in each bin.
    my %binList;
    # Loop through the contigs, processing them.
    for my $contig (@contigs) {
        $stats->Add(contigIn => 1);
        # Is this contig in a bin?
        my $bin = $binMap{$contig};
        if (! $bin) {
            # No, just count it.
            $stats->Add(contigsNotInBins => 1);
        } else {
            # Yes, track the bin.
            $bins{$bin}++;
            push @{$binList{$bin}}, $contig;
            $stats->Add(contigTracked => 1);
        }
    }
    my @binsFound = sort { $bins{$b} <=> $bins{$a} } keys %bins;
    my $binCount = scalar @binsFound;
    if ($binCount == 0) {
        $stats->Add(binlessCluster => 1);
    } else {
        # The cluster intersects one or more bins. Display its data.
        print join("\t", $clusterID, $clusterCount, map { "$_\t$bins{$_}" } @binsFound) . "\n";
        $stats->Add("clusterBins-$binCount" => 1);
        # Do we need to output details?
        if ($oh && $binCount >= $details) {
            for my $bin (sort keys %binList) {
                my $contigs = $binList{$bin};
                for my $contig (@$contigs) {
                    print $oh join("\t", $clusterID, $bin, $contig) . "\n";
                }
            }
        }
    }
    # Update the cluster ID.
    $clusterID++;
}
print "\n\n\n" . $stats->Show(); 
