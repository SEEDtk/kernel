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
use File::Copy::Recursive;
use Stats;
use GPUtils;
use GenomeTypeObject;

=head1 Transfer Genome Packages

    package_transfer.pl [ options ] inDir outDir

This script transfers genome packages between directories. Essentially, if a package passes the
specified criteria, its directory will be copied or moved from the old package repository to the
new one.

=head2 Parameters

The positional parameters are the names of source genome package repository and the target genome repository.

The command-line options are the following.

=over 4

=item move

If specified, qualifying packages will be moved instead of copied, deleting them from the source directory.

=item missing

Only copy a package if it is not already in the target directory. If C<move> is also specified, qualifying
packages will be deleted if they are already in the target directory.

=item samples

If specified, must be the name of a file containing sample IDs. Only packages from the specified samples will
qualify. The file should be tab-delimited with the sample IDs in the first column.

=item original

If specified, only packages from samples will qualify.

=item good

If specified, only packages with a fine SciKit score of 85 or more and a CheckM completeness score of 80 or
more and a properly-annotated Phenylalamine tRNA synthetase will qualify. A package without a
C<quality.tbl> file will automatically not qualify.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outDir',
        ['move', 'delete qualifying packages from the input directory'],
        ['missing', 'do not copy packages already in the output directory'],
        ['samples=s', 'list of samples required for qualification'],
        ['good', 'only good packages qualify'],
        ['original', 'only packages from unmodified bins qualify'],
        );
my $stats = Stats->new();
# Get the input and output directories.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Invalid or missing output directory $outDir.";
}
# Get the options.
my $samples;
if ($opt->samples) {
    $samples = {};
    open(my $ih, '<', $opt->samples) || die "Could not open samples file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\S+)/) {
            $samples->{$1} = 1;
        }
    }
}
my $move = $opt->move;
my $missing = $opt->missing;
my $good = $opt->good;
my $original = $opt->original;
# Get the list of incoming packages.
opendir(my $dh, $inDir) || die "Could not open input directory $inDir: $!";
my @inPackages = grep { -s "$inDir/$_/data.tbl" } readdir $dh;
closedir $dh; undef $dh;

# Get a hash of packages already in the output directory.
opendir($dh, $outDir) || die "Could not open output directory $outDir: $!";
my %outPackages = map { $_ => 1 } grep { substr($_,0,1) ne '.' && -d "$outDir/$_" } readdir $dh;
closedir $dh; undef $dh;
# Loop through the input directory.
for my $package (sort @inPackages) {
    $stats->Add(inPackages => 1);
    # Get this package's sample.
    my $sampleName = '';
    open(my $ih, '<', "$inDir/$package/data.tbl") || die "Could not open $package metadata: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^Sample Name\t(.+)/) {
            $sampleName = $1;
        }
    }
    # Get this package's quality.
    my $isGood;
    if (open(my $qh, '<', "$inDir/$package/quality.tbl")) {
        $stats->Add(qualityFileFound => 1);
        my @flds = split /\t/, <$qh>;
        # Note that the fields may be blank if we don't have quality data. This
        # automatically disqualifies the package.
        if ($flds[12] && $flds[13] && $flds[12] >= 85 && $flds[13] >= 80) {
            # Here we are probably good. Read the GTO and check the annotations.
            my $gtoFile = "$inDir/$package/bin.gto";
            if (-s $gtoFile) {
                my $gto = GenomeTypeObject->create_from_file("$inDir/$package/bin.gto");
                my $flist = GPUtils::role_to_features($gto, 'Phenylalanyl-tRNA synthetase alpha chain');
                if (scalar @$flist == 1) {
                    my $aa = $flist->[0]{protein_sequence};
                    my $aaLen = length $aa;
                    if ($aaLen >= 209 && $aaLen <= 405) {
                        $isGood = 1;
                    }
                }
            }
        }
    }
    # Compute the sample's properties.
    my $isOriginal = ($sampleName && lc $sampleName ne 'derived');
    my $isMissing = ! $outPackages{$package};
    my $passesFilter = (! $samples || $samples->{$sampleName});
    $stats->Add(packageOriginal => 1) if $isOriginal;
    $stats->Add(packageMissing => 1) if $isMissing;
    $stats->Add(packageFiltersIn => 1) if $passesFilter;
    $stats->Add(packageGood => 1) if $isGood;
    # Determine if this sample qualifies.
    if (($isOriginal || ! $original) && ($isGood || ! $good) && $passesFilter) {
        # Here the package qualifies.
        $stats->Add(packageQualifies => 1);
        # Figure out what we should do.
        if ($missing && ! $isMissing) {
            # Here we are only transferring directories not in the target, but this package is in
            # the target. If we are moving, delete the source. Otherwise, do nothing.
            if ($move) {
                $stats->Add(sourceDeleted => 1);
                print "Erasing $inDir/$package.\n";
                File::Copy::Recursive::pathrmdir("$inDir/$package");
            } else {
                print "Package already in target-- skipped.\n";
                $stats->Add(sourceNotCopied => 1);
            }
        } else {
            # Here we are going ahead with a transfer.
            if (! $isMissing) {
                # Delete the target package to make room for the source.
                $stats->Add(targetDeleted => 1);
                print "Deleting $outDir/$package.\n";
                File::Copy::Recursive::pathrmdir("$outDir/$package");
            }
            # Either move or copy.
            if ($move) {
                print "Moving $package to $outDir.\n";
                File::Copy::Recursive::dirmove("$inDir/$package", "$outDir/$package");
                $stats->Add(sourceMoved => 1);
            } else {
                print "Copying $package to $outDir.\n";
                File::Copy::Recursive::dircopy("$inDir/$package", "$outDir/$package");
            }
        }
    }
}
print "All done.\n" . $stats->Show();
