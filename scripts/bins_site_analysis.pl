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
use GenomeTypeObject;

=head1 Analyze the Binning Directories by Site

    bins_site_analysis.pl [ options ] binDir

This script will process a binning jobs directory and produce statistics on the binning jobs organized by site. It will
display site name, number of good bins, number of bad bins, total base pairs, binned base pairs, and unbinned base pairs for each
binning job. The output will be sorted by site.

=head2 Parameters

The positional parameter is the name of the directory containing the binning jobs.

The command-line options are the following.

=over 4

=item good

Name of the good-genome directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ['good=s', 'good-genomes directory', { default => "$FIG_Config::data/GoodPackages"}],
        );
# Get the parameters.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No bin job directory specified.";
} elsif (! -d $binDir) {
    die "Bin job directory $binDir missing or invalid.";
}
my $stats = Stats->new();
# Get the good genomes.
print STDERR "Checking good genomes.\n";
opendir(my $dh, $opt->good) || die "Could not open good-genomes directory: $!";
my %goods = map { $_ => 1 } grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
closedir $dh; undef $dh;
print STDERR scalar(keys %goods) . " good genomes found.\n";
# Find the binning directories.
print STDERR "Gathering bin jobs.\n";
opendir($dh, $binDir) || die "Could not open bin job directory $binDir: $!";
my @binJobs = sort grep { -f "$binDir/$_/bins.rast.json" } readdir $dh;
closedir $dh; undef $dh;
print STDERR scalar(@binJobs) . " completed bin jobs found.\n";
# This hash will store a list of binning-data tuples by site. Each tuple will be
# [name, bpTotal, #good, #bad, bpBinned, bpUnbinned, bpCovg].
my %report;
# Loop through the bins.
for my $binJob (@binJobs) {
    $stats->Add(binFound => 1);
    my $jobDir = "$binDir/$binJob";
    # Initialize our stats for this bin.
    my ($site, $good, $bad, $bpTotal, $bpBinned) = ('(Unspecified)', 0, 0, 0, 0);
    if (! -s "$jobDir/site.tbl") {
        $stats->Add(jobSiteless => 1);
    } else {
        # Check the site.
        open(my $sh, "<$jobDir/site.tbl") || die "Could not open site file: $!";
        my $line = <$sh>;
        (undef, $site) = split /\t/, $line;
    }
    print STDERR "Processing $binJob.\n";
    # This will be set to TRUE when we have processed all the bins.
    my $done;
    # This will be the current bin number.
    my $idx = 1;
    # Loop through the bins.
    while (! $done) {
        # Check for a bin GTO.
        if (-s "$jobDir/bin$idx.gto") {
            # We have one. Extract the ID and check for quality.
            my $gto = GenomeTypeObject->create_from_file("$jobDir/bin$idx.gto");
            my $id = $gto->{id};
            if ($goods{$id}) {
                $good++;
                $stats->Add(binsGood => 1);
            } else {
                $bad++;
                $stats->Add(binsBad => 1);
            }
            # Count the contigs.
            for my $contig (@{$gto->{contigs}}) {
                $bpBinned += length($contig->{dna});
            }
            # Get the next one.
            $idx++;
        } else {
            # We don't have one, so we are done.
            $done = 1;
        }
    }
    if ($good + $bad == 0) {
        $stats->Add(jobBinless => 1);
    } else {
        $stats->Add(jobWithBins => 1);
    }
    # Now count the total base pairs for the sample and compute the coverage and average length.
    my ($covg, $len, $count) = (0, 0, 0);
    open(my $ih, "<$jobDir/contigs.fasta") || die "Could not open contigs.fasta for $jobDir: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        if ($line =~/^>NODE_\d+_length_(\d+)_cov_(\d+\.\d+)/) {
            $stats->Add(spadesContigsRead => 1);
            my ($len0, $covg0) = ($1, $2);
            $covg += $covg0 * $len0;
            $len += $len0;
            $count++;
        } elsif ($line =~ /^>/) {
            $stats->Add(sampleContigsRead => 1);
        } else {
            $bpTotal += length($line);
        }
    }
    # Compute the coverage and mean length.
    my ($avgCovg, $avgLen) = ("n/k", "n/k");
    if ($count > 0 && $len > 0) {
        $avgCovg = $covg / $len;
        $avgLen = $len / $count;
    }
    # Save the results.
    push @{$report{$site}}, [$binJob, $bpTotal, $good, $bad, $bpBinned, $bpTotal - $bpBinned, $avgCovg, $avgLen];
}
# Now we process the output, one site at a time.
print STDERR "Producing reports.\n";
print join("\t", qw(site name size good bad binned unbinned coverage) ) . "\n";
for my $site (sort keys %report) {
    # Sort the site data by total base pairs.
    my @tuples = sort { $a->[1] <=> $b->[1] } @{$report{$site}};
    for my $tuple (@tuples) {
        # Print this job.
        print join("\t", $site, @$tuple) . "\n";
    }
}
print STDERR "All done.\n" . $stats->Show();
