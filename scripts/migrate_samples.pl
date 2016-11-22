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

=head1 Migrate Samples to the Processing Directory

    migrate_samples.pl [ options ] inDir outDir

This script copies samples to the processing directory so they can be submitted to the pipeline. It is a
dirt-simple script that simply copies the first N directories not already present.

=head2 Parameters

The positional parameters are the name of the staging directory for samples followed by the name of the processing
directory.

The command-line options are the following.

=over 4

=item max

Maximum number of directories to copy. The default is C<6>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outDir',
        ['max=i', 'maximum number of samples to move', { default => 6 }],
        );
# Get the directories.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "Input directory missing.";
} elsif (! -d $inDir) {
    die "Invalid input directory $inDir.";
} elsif (! $outDir) {
    die "Output directory missing.";
} elsif (! -d $outDir) {
    die "Invalid output directory $outDir.";
}
# Get the input samples.
my $dh;
opendir($dh, $inDir) || die "Could not open input directory: $!";
my @samples = sort grep { substr($_,0,1) ne '.' && -d "$inDir/$_" } readdir $dh;
closedir $dh; undef $dh;
print scalar(@samples) . " input samples found.\n";
# Get the existing samples.
opendir($dh, $outDir) || die "Could not open output directory: $!";
my %existing = map { $_ => 1 } grep { substr($_,0,1) ne '.' && -d "$outDir/$_" } readdir $dh;
closedir $dh; undef $dh;
print scalar(keys %existing) . " samples already in output directory.\n";
my $moved = 0;
my $max = $opt->max;
print "Processing input samples.\n";
for my $sample (@samples) { last if $moved >= $max;
    if ($existing{$sample}) {
        print "Skipping $sample.\n";
    } else {
        print "Moving $sample to $outDir. ";
        my $numCopied = File::Copy::Recursive::dircopy("$inDir/$sample", "$outDir/$sample");
        if ($numCopied) {
            print "$numCopied items transferred.\n";
            $moved++;
        } else {
            die "Error in copy: $!";
        }
    }
}
print "All done. $moved directories moved.\n";