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

=head1 Create Table of Bin Data

    bins_table.pl [ options ] directory

This program generates a table of bins from a specified directory. Only completed binning projects
will be listed.

A binning project is considered completed if it has a C<bins.rast.json> file. For each such project,
we will iterate through the bins (which are always numbered sequentially from C<1>), and list (0) the
sample ID (which is the directory name), (1) the bin number, (2) the pipeline name (bins_generate),
(3) the name and path of the contig file, and (4) the name and path of the GTO file.

=head2 Parameters

There are two positional parameters-- the name of the directory containing the binning projects, and the
name to give to the output file. Each project is a subdirectory of the specified binning directory.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory outfile');
my $stats = Stats->new();
# Get the main directory name.
my ($directory, $outFile) = @ARGV;
if (! $directory) {
    die "No directory specified.";
} elsif (! -d $directory) {
    die "Invalid directory name $directory.";
}
# Open the output file.
open(my $oh, ">$outFile") || die "Could not open output file: $!";
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    if (! -s "$subDir/bins.rast.json") {
        print "$dir not completed-- skipped.\n";
        $stats->Add(dirsSkipped => 1);
    } else {
        $stats->Add(dirsKept => 1);
        # Now we loop through the bins.
        my $bin = 1;
        while (-s "$subDir/bin$bin.fa") {
            $stats->Add(binsFound => 1);
            print join("\t", $dir, $bin, 'bins_generate', "$subDir/bin$bin.fa", "$subDir/bin$bin.gto") . "\n";
            $bin++;
        }
    }
}
print "All done:\n" . $stats->Show();
