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

=head1 Find the Binning Directories for a Given Site

    bins_site_analysis.pl [ options ] site_name binDir

This script will search a binning jobs directory for bins from a given site with one or more genome results.
For each bin, the sample ID and the number of good and bad output genomes will be displayed.

=head2 Parameters

The positional parameters are the name of the site of interest and the directory containing the binning jobs.

The command-line options are the following.

=over 4

=item good

Name of the good-genome directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('site_name binDir',
        ['good=s', 'good-genomes directory', { default => "$FIG_Config::data/GoodPackages"}],
        );
# Get the parameters.
my ($siteName, $binDir) = @ARGV;
if (! $siteName) {
    die "No site name specified.";
} elsif (! $binDir) {
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
print STDERR scalar(@binJobs) . " bin jobs found.\n";
# Print a header.
print "Sample\tGood\tBad\n";
# Loop through the bins.
for my $binJob (@binJobs) {
    $stats->Add(binFound => 1);
    my $jobDir = "$binDir/$binJob";
    if (! -s "$jobDir/site.tbl") {
        $stats->Add(binSiteless => 1);
    } else {
        # Check the site.
        open(my $sh, "<$jobDir/site.tbl") || die "Could not open site file: $!";
        my $line = <$sh>;
        my ($proj, $site) = split /\t/, $line;
        if ($site ne $siteName) {
            # The bin is for the wrong site.
            $stats->Add(jobWrongSite => 1);
        } elsif (-z "$jobDir/bins.rast.json") {
            # Here the bin job is for the site, but has no bins.
            $stats->Add(jobNoGenomes => 1);
        } else {
            # Here we have a bin job in the correct site which produced genomes.
            $stats->Add(jobWithGenomes => 1);
            print STDERR "Processing $binJob.\n";
            # This will be set to TRUE when we have processed all the bins.
            my $done;
            # This will be the current bin number.
            my $idx = 1;
            # This will count the good and bad bins.
            my ($good, $bad);
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
                    # Get the next one.
                    $idx++;
                } else {
                    # We don't have one, so we are done.
                    $done = 1;
                }
            }
            print "$binJob\t$good\t$bad\n";
        }
    }
}
print STDERR "All done.\n" . $stats->Show();
