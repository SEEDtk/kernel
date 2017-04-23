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
use GPUtils;
use File::Copy::Recursive;
use GenomeTypeObject;

=head1 Bad Package Remover

    package_fixup.pl [ options ] packageDir sampleDir

This script searches for damaged packages caused by the bad-ref-genome bug. We are looking for packages that are
(1) missing a C<bin.gto> file or (2) have a C<bin.gto> file that does not match the C<bin>I<N>C<.gto> file in the
Any bad packages found will be deleted.

=head2 Parameters

The positional parameters are the name of the package directory and the name of the sample directory.

The command-line options are the following.

=over 4

=item test

If specified, bad packages will be identified but not deleted.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir sampleDir',
        ['test', 'do not delete bad packages']
        );
# Get the input directories.
my ($packageDir, $sampleDir) = @ARGV;
if (! $packageDir) {
    die "No package directory specified: $!";
} elsif (! -d $packageDir) {
    die "Missing or invalid package directory $packageDir.";
} elsif (! $sampleDir) {
    die "No sample directory specified: $!";
} elsif (! -d $sampleDir) {
    die "Missing or invalid sample directory $sampleDir.";
}
# Create a statistics object.
my $stats = Stats->new();
# Get all the packages.
my $gHash = GPUtils::get_all($packageDir);
# Loop through them.
print "Processing genomes.\n";
for my $genome (sort keys %$gHash) {
    # Get this package and its directory.
    $stats->Add(packages => 1);
    my $genomeDir = $gHash->{$genome};
    print "Processing $genome in $genomeDir.\n";
    # Check for a gto.
    my $bad;
    if (! -s "$genomeDir/bin.gto") {
        print "$genome has no GTO file.\n";
        $stats->Add(packageNoGto => 1);
        $bad = 1;
    } else {
        # Read the data file.
        my $dHash = GPUtils::get_data($gHash, $genome);
        my ($sample, $bin) = map { $dHash->{$_} } ('Sample Name', 'Bin Number');
        if (! $sample || ! $bin) {
            # No sample info, this means it is not a sample-sourced bin.
            $stats->Add(packageNotSample => 1);
        } elsif (! -d "$sampleDir/$sample") {
            # Sample does not exist.
            $stats->Add(packageMissingSample => 1);
            print "$genome sample $sample is missing.\n";
            $bad = 1;
        } else {
            # Get the GTO file for the sample bin.
            my $sampleGtoFile = "$sampleDir/$sample/bin$bin.gto";
            my $sampleGto = GenomeTypeObject->create_from_file($sampleGtoFile);
            if (! $sampleGto) {
                print "$sampleGtoFile for $genome is malformed or missing.\n";
                $bad = 1;
                $stats->Add(packageBadBinGto => 1);
            } else {
                my $gto = GPUtils::gto_of($gHash, $genome);
                if ($gto->{id} ne $sampleGto->{id}) {
                    print "$sampleGtoFile has different ID from $genome GTO.\n";
                    $bad = 1;
                    $stats->Add(packageWrongBinGto => 1);
                }
            }
        }
        if ($opt->test) {
            print "$genome should be deleted.\n";
        } else {
            print "Deleting $genomeDir.\n";
            File::Copy::Recursive::pathrmdir($genomeDir);
            if (-d $genomeDir) {
                rmdir $genomeDir;
            }
            $stats->Add(packageDeleted => 1);
        }
    }
}
print "All done.\n" . $stats->Show();
