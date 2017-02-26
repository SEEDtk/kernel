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
use RoleParse;
use GenomeTypeObject;
use Stats;
use SeedUtils;

=head1 Remove Bad Contigs from Genome Packages

    gto_remove_bad_contigs.pl [ options ] packageDir outDir

This script runs through genome packages and creates new FASTA files with bad contigs removed. A contig is considered bad
if it does not contain any expected roles from the SciKit quality assessment.

For each incoming genome package, we need access to the L<GenomeTypeObject>, which will be in the C<bin.gto> file,
the list of expected and found roles, which will be in the C<EvalBySciKit/evaluate.out> file, and the package
specifications, which will be in the C<data.tbl> file. The list of good roles will be computed from the C<evaluate.out>
file-- essentially, all those for whom the expected and actual numbers are equal. Each feature in the GTO will be examined,
and if its assignment contains a good role, the contig ID will be memorized. A second pass is then made through the
GTO to create an output FASTA file of only good contigs. The FASTA file will be in the output directory having the same
name as the input genome package directory name with an extension of C<.fa>.

=head2 Parameters

The positional parameters are the name of the input genome package directory and the output directory for the FASTA files.

The command-line options are the following.

=over 4

=item samples

The name of a file containing samples to process. The file must be tab-delimited with the sample names in the first
column. If this parameter is not specified, all genome packages binned from samples will be processed. Otherwise, only
packages from the specified samples will be processed.

=item roles

The name of a tab-delimited file containing role ID mappings. The role ID is in the first column and the checksum for
matching roles in the second column.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir outDir',
        ['samples=s', 'name of a file listing the samples to process'],
        ['roles=s', 'name of role mapping file', { default => "$FIG_Config::global/roles.in.subsystems"}],
        );
# Check the parameters.
my ($packageDir, $outDir) = @ARGV;
if (! $packageDir) {
    die "No input directory specified.";
} elsif (! -d $packageDir) {
    die "Input directory $packageDir not found or invalid.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Output directory not found or invalid.";
}
# Create the statistics object.
my $stats = Stats->new();
# Read in the role map.
open(my $ih, '<', $opt->roles) || die "Could not open role file: $!";
my %roleMap;
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(roleMapIn => 1);
    if ($line =~ /^(\S+)\t(\S+)/) {
        $roleMap{$1} = $2;
        $stats->Add(roleMapStored => 1);
    }
}
close $ih; undef $ih;
# Read in the samples list, if any.
my $samplesH;
if ($opt->samples) {
    open($ih, '<', $opt->samples) || die "Could not open samples file: $!";
    $samplesH = {};
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(sampleIdIn => 1);
        if ($line =~ /^([^\t]+)/) {
            $samplesH->{$1} = 1;
            $stats->Add(sampleIdStored => 1)
        }
    }
    close $ih; undef $ih;
}
# Get the packages with SciKit evaluations in the input directory.
opendir(my $dh, $packageDir) || die "Could not open input directory: $!";
my @packages = grep { -s "$packageDir/$_/EvalBySciKit/evaluate.out" } readdir $dh;
closedir $dh;
# Loop through the packages, filtering for sample bins.
for my $package (@packages) {
    # Compute the full package directory name.
    my $thisDir = "$packageDir/$package";
    # Check the sample name.
    open(my $fh, '<', "$thisDir/data.tbl") || die "Could not open data.tbl file for 4package: $!";
    $stats->Add(packageChecked => 1);
    my $sampleID;
    while (! eof $fh) {
        if (<$fh> =~ /^Sample Name\t(.+)/) {
            $sampleID = $1;
        }
    }
    if (! $sampleID) {
        $stats->Add(packageNotSample => 1);
    } elsif ($sampleID eq 'Derived') {
        $stats->Add(packageDerived => 1);
    } elsif ($samplesH && ! $samplesH->{$sampleID}) {
        $stats->Add(packageFilteredOut => 1);
    } else {
        print "Processing $package.\n";
        # This hash will contain checksums for all the roles we consider valid.
        my %goodRoleChecksums;
        # Get the SciKit evaluation data for this package.
        open(my $ih, '<', "$thisDir/EvalBySciKit/evaluate.out") || die "Could not open evaluation output for $package: $!";
        # Extract the good roles and save their checksums.
        while (! eof $ih) {
            if (<$ih> =~ /^(\S+)\t(\d+)(?:\.\d+)?\t(\d+)/) {
                my ($roleID, $expect, $actual) = ($1, $2, $3);
                $stats->Add(EvaluateRoleIn => 1);
                if ($expect > 0 && $expect == $actual) {
                    my $checksum = $roleMap{$roleID};
                    if (! $checksum) {
                        print "Role $roleID not found in role file.\n";
                        $stats->Add(missingRole => 1);
                    } else {
                        $stats->Add(EvaluateRoleKept => 1);
                        $goodRoleChecksums{$checksum} = 1;
                    }
                }
            }
        }
        # Read in the GTO.
        my $gto = GenomeTypeObject->create_from_file("$thisDir/bin.gto");
        # This hash will remember the contigs with good roles in them.
        my %goodContigs;
        # Loop through the features.
        print "Processing features.\n";
        my $features = $gto->{features};
        for my $feature (@$features) {
            my $function = $feature->{function};
            if ($function) {
                $stats->Add(functionIn => 1);
                my @roles = SeedUtils::roles_of_function($function);
                for my $role (@roles) {
                    $stats->Add(roleIn => 1);
                    my $checksum = RoleParse::Checksum($role);
                    if ($goodRoleChecksums{$checksum}) {
                        # We have a good role. Get its contigs.
                        $stats->Add(goodRole => 1);
                        my $location = $feature->{location};
                        for my $region (@$location) {
                            $goodContigs{$region->[0]}++;
                            $stats->Add(goodContigFlagged => 1);
                        }
                    }
                }
            }
        }
        # Create the FASTA output.
        my $outFile = "$outDir/$package.fa";
        print "Generating $outFile.\n";
        open(my $oh, '>', $outFile) || die "Could not open FASTA output $outFile: $!";
        my $contigs = $gto->{contigs};
        for my $contig (@$contigs) {
            $stats->Add(contigsIn => 1);
            my $contigID = $contig->{id};
            my $dna = $contig->{dna};
            my $dnaLen = length $dna;
            if ($goodContigs{$contigID}) {
                $stats->Add(contigsOut => 1);
                $stats->Add(dnaOut => $dnaLen);
                print $oh ">$contigID\n$contig->{dna}\n";
            } else {
                $stats->Add(contigsRemoved => 1);
                $stats->Add(dnaRemoved => $dnaLen);
            }
        }
        close $oh;
    }
}
