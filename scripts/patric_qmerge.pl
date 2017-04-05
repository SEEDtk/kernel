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

=head1 Merge PATRIC Quality Data

    patric_qmerge.pl [ options ] qfile1 qfile2 ...

This script reads SciKit quality data from one or more files and merges it into a standard quality report file
(such as produced by L<package_report.pl>). The quality report file has a sample ID in the first column, which
must by C<PATRIC>. The second column is the genome ID, the third is the genome name, and the coarse and fine
SciKit scores go in columns 12 and 13, respectively. The merged quality report file will appear as the standard
output. Status messages will appear on the standard error output.

=head2 Parameters

The positional parameters are the names of the files containing the SciKit quality data. These files must be tab-delimited,
with each record containing (0) a genome ID, (1) a coarse SciKit score, and (2) a fine SciKit score.

The command-line options are those found in L<ScriptUtils/ih_options> (modifying the standard input) plus the following

=over 4

=item lostGenomes

If specified, the name of a file to contain the IDs and names of genomes from the input file for which no valid SciKit
score could be found. If scores are already present in the input file they will not be included in this list.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('qfile1 qfile2 ...', ScriptUtils::ih_options(),
        ['lostGenomes=s', 'name of file for genomes whose scores are needed']
        );
my $stats = Stats->new();
# We need a hash of the available SciKit scores.
my %scores;
# Loop through the quality files.
my @files = @ARGV;
print STDERR scalar(@files) . " quality files specified.\n";
for my $qFile (@files) {
    open(my $qh, "<$qFile") || die "Could not open $qFile: $!";
    $stats->Add(qFileIn => 1);
    print STDERR "Scanning $qFile.\n";
    while (! eof $qh) {
        my $line = <$qh>;
        $stats->Add(qFileLineIn => 1);
        if ($line =~ /^(\d+\.\d+)\t(\d+.\d+)\t(\d+.\d+)$/) {
            my ($genomeID, $skC, $skF) = ($1, $2, $3);
            if ($skC > 100 || $skF > 100) {
                $stats->Add(qFileBadScores => 1);
            } else {
                $scores{$genomeID} = [$skC, $skF];
                $stats->Add(qFileScore => 1);
            }
        } else {
            $stats->Add(qFileBadLine => 1);
        }
    }
}
print STDERR "Scores found for " . scalar(keys %scores) . " genomes.\n";
# Check for a lost-genomes option.
my $fh;
if ($opt->lostgenomes) {
    open($fh, '>', $opt->lostgenomes) || die "Could not open lostGenomes file: $!";
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
print STDERR "Processing input file.\n";
# Loop through the genomes.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    $stats->Add(inputLineIn => 1);
    # Parse the input line.
    my @fields = split /\t/, $line;
    my ($source, $genome, $name) = @fields;
    if ($source eq 'PATRIC') {
        $stats->Add(patricLineIn => 1);
        if (! $scores{$genome}) {
            # Here we do not have scores to put in.
            $stats->Add(sciKitNotFound => 1);
            if ($fields[11] == 0 && $fields[12] == 0) {
                # Here we need scores and do not have them.
                $stats->Add(lostGenome => 1);
                if ($fh) {
                    print $fh "$genome\t$name\n";
                }
            }
        } else {
            # Here we have new scores.
            $stats->Add(sciKitFound => 1);
            ($fields[11], $fields[12]) = @{$scores{$genome}};
        }
    }
    # Output the data line.
    print join("\t", @fields) . "\n";
    $stats->Add(lineOut => 1);
}
print STDERR "All done.\n" . $stats->Show();