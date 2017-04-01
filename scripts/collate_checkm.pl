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

=head1 Organize CheckM output files.

    collate_checkm.pl [ options ] inDir1 inDir2 ... inDirN

This script reads one or more directories for files with a suffix of C<.out>. Such files
are presumed to contain CheckM output. The output will be organized into a single tab-delimited
output file with the columns (0) bin ID, (1) marker taxonomy, (2) completeness, (3) contamination,
and (4) strain heterogeneity. CheckM output files use multiple space characters as column delimiters,
with delimiters at each end. The columns of interest are in position 0, 1, 11, 12, and 13.

=head2 Parameters

The positional parameters are the names of the directories to process.

The following command-line options are supported.

=over 4

=item altFile

If specified, the name of a file containing alternative quality data. The file should be tab-delimited, with
the genome ID in the first column. The CheckM values will be added at the end.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir1 inDir2 ... inDirN',
        ["altFile=s", 'alternative quality data file'],
        );
# Get the directory names.
my @inDirs = @ARGV;
# Check for an alternative quality file.
my %qual;
if ($opt->altfile) {
    open(my $qh, "<", $opt->altfile) || die "Could not open altFile: $!";
    while (! eof $qh) {
        my $line = <$qh>;
        if ($line =~ /^(\d+\.\d+)\t(.+)/) {
            $qual{$1} = $2;
        }
    }
}
# Loop through the directories.
for my $inDir (@inDirs) {
    if (! -d $inDir) {
        die "$inDir is not a valid directory.";
    } elsif (! opendir(my $dh, $inDir)) {
        die "Could not open $inDir: $!";
    } else {
        my @logFiles = grep { substr($_, -4) eq '.out' } readdir $dh;
        # Loop through the files in the directory.
        for my $logFile (@logFiles) {
            open(my $ih, "<$inDir/$logFile") || die "Could not open $logFile: $!";
            # Loop through the lines of the file. We look for lines beginning with genome IDs.
            while (! eof $ih) {
                my $line = <$ih>;
                if ($line =~ /^\s\s+(\d+\.\d+)\s\s+(.+)/) {
                    my ($genome, $data) = ($1, $2);
                    my @fields = split /\s\s+/, $data;
                    my ($lineage, undef, undef, undef, undef, undef, undef,
                        undef, undef, undef, $completeness, $contamination, $heterogeneity) =
                        @fields;
                    my $qual = $qual{$genome};
                    if (! $opt->altfile || $qual) {
                        if ($qual) {
                            $genome = "$genome\t$qual";
                        }
                        print join("\t", $genome, $lineage, $completeness, $contamination, $heterogeneity) . "\n";
                    }
                }
            }
        }
    }
}
