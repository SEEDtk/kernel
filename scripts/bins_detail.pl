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
use Shrub;
use Bin;
use Bin::Analyze;
use File::Copy::Recursive;
use Loader;

=head1 Evaluate Community Bins

    bins_detail.pl [ options ] binDir

This program displays a report on the contigs in each bin, listing the closest genome, tetranucleotide distance, and coverage
distance for each genome in the bin.

=head2 Parameters

There is one positional parameter-- the directory containing the bin information. The json file containing the bin data structures
must be named C<bins.json> and the bin data for the contigs must be in C<contigs.ref.bins>. These are output by
L<bins_create.pl>.

The command-line options are those found in L<Shrub/script_options>.


=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDirectory', Shrub::script_options(),
        );
# Create the loader object and get the statistics.
my $loader = Loader->new();
my $stats = $loader->stats;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Read the bins.
my $binList;
my ($binDir) = @ARGV;
if (! $binDir) {
    die "A bins directory is required.";
} elsif (! -d $binDir) {
    die "Bin directory $binDir not found or invalid.";
}
# Compute the bins file name.
my $jsonFileName = "$binDir/bins.json";
print "Reading bins.\n";
$binList = Bin::ReadBins($jsonFileName);
# Read in the single-contig bin file.
my $cbinList = Bin::ReadContigs("$binDir/contigs.ref.bins");
my %contigBins = map { $_->contig1 => $_ } @$cbinList;
# Loop through the bins, producing the report.
my $analyzer = Bin::Analyze->new($shrub);
# Loop through the bins.
for (my $id = 0; $id < scalar(@$binList); $id++) {
    my $bin = $binList->[$id];
    $stats->Add(binsProcessed => 1);
    print "Processing bin $id.\n";
    $analyzer->ContigReport(\*STDOUT, $id, $bin, \%contigBins);
}
print "All done.\n" . $stats->Show();