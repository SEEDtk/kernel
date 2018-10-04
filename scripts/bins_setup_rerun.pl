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

=head1 Set Up to Rerun Binning

    bins_setup_rerun.pl [ options ] oldDir newDir

This is a simple script that transfers binning jobs to a new directory so they can be rerun. The binning directories are created
in the new location with the same names they previously had. Only the key input files are transferred-- C<contigs.fasta>,
C<output.contigs2reads.txt>, C<site.tbl>, and I<name>C<_abundance_table.tsv>.  An empty C<Assembly> subdirectory is created
so that the binning job looks assembled but not cleaned.

=head2 Parameters

The positional parameters are the name of the source directory and the name of the new directory. If the new directory does not
exist, it will be created.

The following command-line options are supported.

=over 4

=item clear

Erase the output directory before starting. If omitted, binning jobs that already exist in the output will not be copied.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('oldDir newDir',
        ['clear', 'erase new directory before starting']);
# Get the specified directories.
my ($oldDir, $newDir) = @ARGV;
if (! $oldDir) {
    die "No input directory specified.";
} elsif (! -d $oldDir) {
    die "Input directory $oldDir not found.";
} elsif (! $newDir) {
    die "No output directory specified.";
} elsif (! -d $newDir) {
    print "Creating $newDir.\n";
    File::Copy::Recursive::pathmk($newDir) || die "Could not create $newDir: $!";
} elsif ($opt->clear) {
    print "Erasing $newDir.\n";
    File::Copy::Recursive::pathempty($newDir) || die "Could not clear $newDir: $!";
}
my $stats = Stats->new();
# Read the input directory.
print "Scanning $oldDir.\n";
opendir(my $dh, $oldDir) || die "Could not open $oldDir: $!";
my @jobs = grep { -s "$oldDir/$_/contigs.fasta" } readdir $dh;
close $dh;
my $total = scalar(@jobs);
print "$total binning jobs found in $oldDir.\n";
my $count = 0;
# Loop through the binning jobs.
for my $job (@jobs) {
    my $source = "$oldDir/$job";
    my $target = "$newDir/$job";
    $count++;
    if (-d $target) {
        print "$job already in $newDir-- skipped.\n";
        $stats->Add(dirSkipped => 1);
    } else {
        print "Processing $job ($count of $total).\n";
        File::Copy::Recursive::pathmk($target) || die "Could not create $target: $!";
        # Start a list of the files to copy.
        my @files = qw(contigs.fasta output.contigs2reads.txt);
        my $abundance = $job . '_abundance_table.tsv';
        if (-s "$source/$abundance") {
            push @files, $abundance;
            $stats->Add(abundanceFound => 1);
        }
        if (-s "$source/site.tbl") {
            push @files, 'site.tbl';
            $stats->Add(siteFound => 1);
        }
        # Copy the files.
        for my $file (@files) {
            File::Copy::Recursive::fcopy("$source/$file", "$target/$file") || die "Could not copy $file: $!";
            $stats->Add(fileCopied => 1);
        }
        # Create the assembly directory.
        File::Copy::Recursive::pathmk("$target/Assembly");
        $stats->Add(dirCopied => 1);
    }
}
print "All done.\n" . $stats->Show();