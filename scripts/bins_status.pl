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
use File::Copy::Recursive;

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

=item Done

Expectation processing completed-- C<expect.report.txt> exists or C<bins.rast.json> exists and there is no
abundance file.

=back

=head2 Parameters

The single positional parameter is the name of the directory containing the sample sub-directories.

The following command-line options are supported.

=over 4

=item clean

Remove assembly information for completed jobs.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory',
                ['clean', 'clean up assembly information for complete samples']);
my $stats = Stats->new();
# Get the main directory name.
my ($directory) = @ARGV;
if (! $directory) {
    die "No directory specified.";
} elsif (! -d $directory) {
    die "Invalid directory name $directory.";
}
# Save the options.
my $clean = $opt->clean;
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    my $rastFound = (-s "$subDir/bins.rast.json");
    my $cleaned = (-d "$subDir/Assembly" ? "" : "  Assembly cleaned.");
    my $done;
    # Determine the site.
    my $site;
    if (! -s "$subDir/site.tbl") {
        $site = "Unspecified";
    } elsif (open(my $sh, '<', "$subDir/site.tbl")) {
        my $line = <$sh>;
        if ($line =~ /(\S+)\t([^\t]+)\t(.+)/) {
            $site = "$1 $3";
            $stats->Add("site-$2" => 1);
        } else {
            $site = "Invalid";
        }
    } else {
        $site = "Error";
    }
    my $label = "$subDir ($site)";
    # Determine the status.
    if (-s "$subDir/expect.report.txt") {
        $done = "Expectations Computed.";
    } elsif ($rastFound && ! -s "$subDir/$dir" . '_abundance_table.tsv') {
        $done = "Done (No Expectations).";
    } elsif ($rastFound) {
        print "$label: RAST Complete.\n";
        $stats->Add(dirs5RastComplete => 1);
    } elsif (-s "$subDir/bin1.gto") {
        print "$label: RAST in Progress.\n";
        $stats->Add(dirs4RastPartial => 1);
    } elsif (-s "$subDir/bins.json") {
        print "$label: Bins Computed.\n";
        $stats->Add(dirs3Binned => 1);
    } elsif (-s "$subDir/Assembly/contigs.fasta") {
        print "$label: Assembled.\n";
        $stats->Add(dirs2Assembled => 1);
    } elsif (-d "$subDir/Assembly") {
        print "$label: Assembling.\n";
        $stats->Add(dirs1Assembling => 1);
    } else {
        print "$label: Downloaded.\n";
        $stats->Add(dirs0Downloaded => 1);
    }
    # If we are done, we process here and check for cleaning.
    if ($done) {
        $stats->Add(dirs6Done => 1);
        if ($clean && ! $cleaned) {
            $cleaned = "  Cleaning Assembly.";
            File::Copy::Recursive::pathempty("$subDir/Assembly");
            rmdir "$subDir/Assembly";
            opendir(my $dh, $subDir) || die "Could not open $subDir: $!";
            my @fastqs = grep { $_ =~ /\.fastq$/ } readdir $dh;
            for my $fastq (@fastqs) {
                unlink "$subDir/$fastq";
            }
            $stats->Add(dirsCleaned => 1);
        }
        print "$label: $done$cleaned\n";
    }
}
print "\nAll done:\n" . $stats->Show();
