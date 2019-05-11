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

=head1 Audit Community Bins

    bins_audit.pl [ options ] binDir

This program performs an analysis of the progress of the current binning effort.  For each type of sample,
we will accumulate the number of good bins, the number of bad bins, the number of bins in progress,  the number of
failed binning attempts, the number of unprocessed samples, the number of completed samples, and the number of
samples in progress.  The hope is that this information can be used to compute the total number of expected good
bins.

=head2 Parameters

The sole positional parameter is the name of the binning directory.  This must contain samples as subdirectories.

The command-line options are the following.

=over 4

=item verbose

Write progress messages to STDERR.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ['verbose|debug|v', 'write progress to STDERR']
        );
# Create the statistics object.
my $stats = Stats->new();
# Get the binning directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir is invalid or missing.";
}
# Get the options.
my $debug = $opt->verbose;
# This hash will track the totals for each site type.  The columns are
# [site name, good-bins, bad-bins, accum-bins, fails, done, in-progress, not-started].
use constant { SITE_NAME => 0, GOOD_BINS => 1, BAD_BINS => 2, ACCUM_BINS => 3, FAILS => 4,
        DONE => 5, IN_PROGRESS => 6, NOT_STARTED => 7 };
my %cats;
# Loop through the binning directories.
print STDERR "Reading $binDir.\n" if $debug;
opendir(my $dh, $binDir) || die "Could not open binning directory: $!";
while (my $sample = readdir $dh) {
    my $subDir = "$binDir/$sample";
    my $siteFile = "$subDir/site.tbl";
    if (substr($sample, 0, 1) ne '.' && -s $siteFile) {
        $stats->Add(sampleFound => 1);
        open(my $ih, '<', $siteFile) || die "Could not open $siteFile: $!";
        my $line = <$ih>;
        close $ih; undef $ih;
        chomp $line;
        my (undef, undef, $cat) = split /\t/, $line;
        print STDERR "Processing $subDir: $cat\n" if $debug;
        if (! $cats{$cat}) {
            $stats->Add(catFound => 1);
            $cats{$cat} = [$cat, 0, 0, 0, 0, 0, 0, 0];
        }
        my $counters = $cats{$cat};
        # Are we done, in-progress, or not-started?
        if (-s "$subDir/Eval/index.tbl") {
            # Here we are done.  Did we fail to find any bins?
            if (! -s "$subDir/ref.genomes.scores.tbl") {
                $stats->Add(sampleFail => 1);
                $counters->[FAILS]++;
            } else {
                $stats->Add(sampleDone => 1);
                $counters->[DONE]++;
                # Bins were found.  We count them.
                open($ih, '<', "$subDir/Eval/index.tbl") || die "Could not open evaluation file for $sample: $!";
                my $line = <$ih>;
                while (! eof $ih) {
                    $line = <$ih>;
                    if ($line =~ /\t1$/) {
                        $counters->[GOOD_BINS]++;
                        $stats->Add(binsGood => 1);
                    } else {
                        $counters->[BAD_BINS]++;
                        $stats->Add(binsBad => 1);
                    }
                }
            }
        } elsif (-f "$subDir/sample.fasta") {
            # Here we are in-progress.  Count the bins accumulating.
            $stats->Add(sampleInProgress => 1);
            $counters->[IN_PROGRESS]++;
            my $i = 1;
            while (-s "$subDir/bin$i.gto") {
                $i++;
            }
            $i--;
            $stats->Add(binsInProgress => $i);
            $counters->[ACCUM_BINS]++;
        } else {
            # Here we are not started.
            $stats->Add(sampleNotStarted => 1);
            $counters->[NOT_STARTED]++;
        }
    }
}
print STDERR "Producing report.\n" if $debug;
print join("\t", qw(site name good_bins bad_bins accum_bins fail_samples done_samples in_progress not_started)) . "\n";
for my $cat (sort keys %cats) {
    print join("\t", @{$cats{$cat}}) . "\n";
}
print STDERR "All done.\n" . $stats->Show() if $debug;