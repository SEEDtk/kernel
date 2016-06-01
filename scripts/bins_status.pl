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

=head1 Check Status of Bin Pipeline

    bins_status.pl [ options ] directory

Check the subdirectories of a specified directory to determine the progress of the binning runs. The status possibilities
are

=over 4

=item Assembling

The reads are being assembled into contigs-- C<Assembly> exists.

=item Assembled

The reads have been assembled into contigs-- C<Assembly/contigs.fasta> exists.

=item Bins Computed

The contigs have been sorted into bins-- C<bins.json> exists.

=item RAST in Progress

One or more bins have been processed by RAST-- C<bin>I<X>C<.gto> exists for one or more bin numbers I<X>.

=item RAST Complete

RAST processing completed-- C<bins.rast.json> exists.

=item Expectations Computed

Expectation processing completed-- C<expect.report.txt> exists.

=back

=head2 Parameters

The single positional parameter is the name of the directory containing the sample sub-directories. 

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory');
my $stats = Stats->new();
# Get the main directory name.
my ($directory) = @ARGV;
if (! $directory) {
    die "No directory specified.";
} elsif (! -d $directory) {
    die "Invalid directory name $directory.";
}
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    # Determine the status.
    if (-s "$subDir/expect.report.txt") {
        print "$subDir: Expectations Computed.\n";
        $stats->Add(dirs6Expect => 1);
    } elsif (-s "$subDir/bins.rast.json") {
        print "$subDir: RAST Complete.\n";
        $stats->Add(dirs5RastDone => 1);
    } elsif (-s "$subDir/bin1.gto") {
        print "$subDir: RAST in Progress.\n";
        $stats->Add(dirs4RastPartial => 1);
    } elsif (-s "$subDir/bins.json") {
        print "$subDir: Bins Computed.\n";
        $stats->Add(dirs3Binned => 1);
    } elsif (-s "$subDir/Assembly/contigs.fasta") {
        print "$subDir: Assembled.\n";
        $stats->Add(dirs2Assembled => 1);
    } elsif (-d "$subDir/Assembly") {
        print "$subDir: Assembling.\n";
        $stats->Add(dirs1Assembling => 1);
    } else {
        print "$subDir: Downloaded.\n";
        $stats->Add(dirs0Downloaded => 1);
    }
}
print "\nAll done:\n" . $stats->Show();
