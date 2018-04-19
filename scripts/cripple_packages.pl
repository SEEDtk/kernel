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
use GenomeTypeObject;
use File::Copy::Recursive;
use GPUtils;

=head1 Create Packages By Crippling GTOs

    cripple_packages.pl [ options ] inDir outDir

This script removes a certain percentage of features from one or more genome packages, producing new files in
a different directory. It does not alter the DNA sequences in any way, only the annotations. This is used for
quality tool testing.

We remove features randomly, so multiple runs of this script will generate different results.

If a package already exists in the output directory for a specific genome, it will not be overwritten.

=head2 Parameters

The positional parameters are the names of the input package directory and the output package directory.

The command-line options are the following.

=over 4

=item removal

Percent of features to remove. The default is C<25>.

=item clear

If specified, the output directory will be cleared before processing.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outDir',
        ['removal|r=f', 'percent of features to remove', { default => 25 }],
        ['clear', 'clear output directory before starting']
        );
# Compute the removal percentage.
my $removal = $opt->removal;
# Get the parameters.
my ($inDir, $outDir) = @ARGV;
# Check the directories.
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Input directory $inDir not found or invalid.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} elsif ($opt->clear) {
    print "Clearing output directory $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not clear $outDir: $!";
}
# Get the list of input packages.
print "Reading $inDir directory.\n";
my $gHash = GPUtils::get_all($inDir);
# Loop through the packages.
for my $genome (sort keys %$gHash) {
    # Compute the output directory path.
    my $inPath = $gHash->{$genome};
    my $outPath = File::Spec->catfile($outDir, File::Spec->abs2rel($inPath, $inDir));
    if (-d $outPath) {
        print "$genome is already in $outDir-- skipping.\n";
    } else {
        print "Processing $genome.\n";
        File::Copy::Recursive::pathmk($outPath) || die "Could not create new package directory: $!";
        File::Copy::Recursive::fcopy("$inPath/bin.fa", "$outPath/bin.fa") || die "Could not copy bin.fa: $!";
        # Create the new data table.
        my $dataHash = GPUtils::get_data($gHash, $genome);
        $dataHash->{'Cripple Percent'} = $removal;
        # Write the data file.
        open(my $oh, '>', "$outPath/data.tbl") || die "Could not open data table file: $!";
        for my $key (sort keys %$dataHash) {
            print $oh "$key\t$dataHash->{$key}\n";
        }
        close $oh;
        # Get the GTO.
        my $gto = GPUtils::gto_of($gHash, $genome);
        # Remove features, then save it to the target.
        $gto->cripple($removal);
        my $gtoFile = File::Spec->catfile($outPath, 'bin.gto');
        print "Writing crippled GTO to $gtoFile.\n";
        $gto->destroy_to_file($gtoFile);
    }
}
print "All done.\n";
