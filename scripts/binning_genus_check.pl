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
use GPUtils;

=head1 Display Reference Genome Distribution by Binning Site

    binning_genus_check.pl [ options ] binDir pkgDir

This script runs through all of the packages in a directory. For each one, it will read the C<data.tbl> file to determine the genus of the bin
and the source package. It will then look up the site of the package in the specified binning directory using the sample's C<site.tbl> file.
The output will be a count of genus occurrences by site.

=head2 Parameters

The two positional parameters are the directory containing the binning job results and the directory containing the genome packages.


=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir pkgDir',
        );
# Get the two directories.
my ($binDir, $pkgDir) = @ARGV;
if (! $binDir) {
    die "No binning job directory specified.";
} elsif (! -d $binDir) {
    die "Invalid binning job directory $binDir.";
} elsif (! $pkgDir) {
    die "No package directory specified.";
} elsif (! -d $pkgDir) {
    die "Invalid package directory $pkgDir.";
}
# This hash will map sample names to sites.
my %sampleSite;
# This is a hash of counts, two-dimensional by genus within site.
my %siteGenus;
# Get all the package directories.
print STDERR "Reading directory $pkgDir.\n";
my $genomeHash = GPUtils::get_all($pkgDir);
# Loop through the genomes.
for my $genome (sort keys %$genomeHash) {
    print STDERR "Processing $genome.\n";
    # Get its data file.
    my $dataHash = GPUtils::get_data($genomeHash, $genome);
    # Isolate the sample and the genus.
    my $sample = $dataHash->{'Sample Name'};
    my $name = $dataHash->{'Genome Name'};
    my ($genus) = split /\s+/, $name;
    # Compute the site for the sample.
    my $site = $sampleSite{$sample};
    if (! $site) {
        # Here we have to find the site.
        my $siteFile = "$binDir/$sample/site.tbl";
        if (! -s $siteFile) {
            $site = "<unspecified>";
        } else {
            open(my $sh, "<$siteFile") || die "Could not open $siteFile: $!";
            my $siteLine = <$sh>;
            (undef, $site) = split /\t/, $siteLine;
        }
        $sampleSite{$sample} = $site;
    }
    # Now we have the site. Update the count.
    $siteGenus{$site}{$genus}++;
}
# All done. Output the counts.
print STDERR "Writing counts.\n";
print join("\t", qw(site genus count)) . "\n";
for my $site (sort keys %siteGenus) {
    my $genusH = $siteGenus{$site};
    for my $genus (sort keys %$genusH) {
        print join("\t", $site, $genus, $genusH->{$genus}) . "\n";
    }
}
print STDERR "All done.\n";
