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
use Bin::Package;

=head1 Create Genome Packages

    package_bins.pl [ options ] sampleDir packageDir

This script will create genome packages from the bins in a single sample directory. Each genome package is a directory with the
same name as the bin's genome ID and contains the following files

=over 4

=item bin.gto

A L<GenomeTypeObject> for the genome formed from the bin.

=item bin.fa

A FASTA file containing the bin's contigs.

=item data.tbl

A tab-delimited file of key/value pairs, with the following keys.

=over 8

=item Genome Name

Name of the bin genome

=item Sample Name

Name of the sample from which the bin came.

=item Bin Number

The ID number of the bin in the sample.

=item Contigs

The number of contigs in the bin.

=item Base pairs

The number of DNA base pairs in the bin's contigs.

=item Ref Genome

The ID of the closest reference genome.

=item Ref Name

The name of the closest reference genome.

=back

=back

=head2 Parameters

The two positional parameters are the full name of the sample directory and the full name of the directory to contain the output
genome packages. If the genome package already exists, it will not be overwritten.

The command-line options are the following.

=over 4

=item force

If specified, existing genome packages will be overwritten.

=item recursive

If specified, then the sample directory is treated as a directory of sample directories, and all samples therein will be processed.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir packageDir',
        ['force', "overwrite existing packages"],
        ['recursive', "process a directory of samples"],
        );
my $stats = Stats->new();
# Get the directories.
my ($sampleDir, $packageDir) = @ARGV;
if (! $sampleDir) {
    die "No sample directory specified.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
if (! $packageDir) {
    die "No package directory specified.";
} elsif (! -d $packageDir) {
    die "Invalid package directory $packageDir";
}
my @samples;
if (! $opt->recursive) {
    push @samples, $sampleDir;
} else {
    opendir(my $dh, $sampleDir) || die "Could not open sample directory: $!";
    push @samples, map { "$sampleDir/$_" } grep { substr($_, 0, 1) ne '.' && -d "$sampleDir/$_" } readdir $dh;
}
for my $sample (@samples) {
    $stats->Add(samples => 1);
    if (! -s "$sample/bins.rast.json") {
        print "Skipping $sample: incomplete.\n";
        $stats->Add(sampleSkipped => 1);
    } elsif (! -s "$sample/bin1.gto") {
        print "Skipping $sample: no bins.\n";
        $stats->Add(sampleEmpty => 1);
    } else {
        # Get the sample name.
        my @pieces = split /[\\\/:]/, $sample;
        my $sampleName = pop @pieces;
        print "Sample name is $sampleName.\n";
        Bin::Package::CreateFromSample($sample, $sampleName, $stats, $opt->force, $packageDir);
    }
}
print "All done.\n" . $stats->Show();


