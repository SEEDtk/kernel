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

=item blacklist

Name of a file containing the names of samples to be excluded.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outDir',
        ['max=i', 'maximum number of samples to move', { default => 6 }],
        ['blacklist=s', 'file of samples to exclude'],
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
# Check for a blacklist.
my %blackList;
if ($opt->blacklist) {
    open(my $ih, '<', $opt->blacklist) || die "Could not open blacklist file: $!";
    while (! eof $ih) {
        my $name = <$ih>;
        chomp $name;
        $blackList{$name} = 1;
    }
}
# Get the input samples.
my $dh;
opendir($dh, $inDir) || die "Could not open input directory: $!";
my @samples = sort grep { substr($_,0,1) eq 'S' && -d "$inDir/$_" } readdir $dh;
closedir $dh; undef $dh;
print scalar(@samples) . " input samples found.\n";
# Get the existing samples.
opendir($dh, $outDir) || die "Could not open output directory: $!";
my %existing = map { $_ => 1 } grep { substr($_,0,1) ne '.' && -d "$outDir/$_" } readdir $dh;
closedir $dh; undef $dh;
print scalar(keys %existing) . " samples already in output directory.\n";
my $moved = 0;
my $remaining = 0;
my $max = $opt->max;
print "Processing input samples.\n";
for my $sample (@samples) {
    if (! $existing{$sample} && ! $blackList{$sample}) {
        if ($moved < $max) {
            my $target = "$outDir/$sample";
            print "Moving $sample...\n";
            my $numCopied = File::Copy::Recursive::dircopy("$inDir/$sample", $target);
            if ($numCopied) {
                print "$numCopied items transferred.\n";
                $moved++;
            } else {
                print "Error in copy: $!\n";
                if (-d $target) {
                    File::Copy::Recursive::pathrmdir($target);
                }
            }
        } else {
            $remaining++;
        }
    }
}
print "All done. $moved directories moved, $remaining remaining.\n";