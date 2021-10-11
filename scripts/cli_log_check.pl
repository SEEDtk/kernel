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

=head1 Create a QzaReport Progress File from the CLI.UTILS Logs 

    cli_log_check.pl [ options ] logDir outFile

This program reads the CLI.UTILS logs listed in the specified directory and creates a progress file
as output.  The basic strategy is to find lines like

	23 representatives kept out of 75 found in sample ERR2730213.

This line gets output as the progress entry

	ERR2730213	23

However, since some samples are processed multiple times, we only want to keep the last version found.

=head2 Parameters

The positional parameters are the name of the input directory and the name of the output progress file.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir outFile',
        );
# Create a statistics object.
my $stats = Stats->new();
my ($inDir, $outFile) = @ARGV;
if (! $inDir) {
	die "No input directory specified.";
} elsif (! -d $inDir) {
	die "$inDir is not a valid directory.";
} elsif (! $outFile) {
	die "No output file specified.";
}
# Get the list of log files.
opendir (my $dh, $inDir) || die "Could not open input directory: $!";
my @files = grep { $_ =~ /^cli\.utils\.[^.]+\.log$/ } readdir $dh;
closedir $dh;
print scalar(@files) . " files found in directory $inDir.\n";
# This hash will hold the samples found.
my %samples;
# Start the output file.
open(my $oh, '>', $outFile) || die "Could not open output file: $!";
print $oh "sample_id\treps\n";
# Loop through the log files, remembering the status line.
for my $file (@files) {
	my $fileName = "$inDir/$file";
	$stats->Add(fileIn => 1);
	print "Processing log file $fileName.\n";
	open(my $ih, '<', $fileName) || die "Could not open log file $fileName: $!";
	while (! eof $ih) {
		$stats->Add(lineIn => 1);
		my $line = <$ih>;
		if ($line =~ /(\d+) representatives kept out of \d+ found in sample (\S+)\./) {
			$samples{$2} = $1;
			$stats->Add(sampleFound => 1);
		}
	}
}
# Unspool the memorized sample IDs.
for my $sample (sort keys %samples) {
	print $oh "$sample\t$samples{$sample}\n";
	$stats->Add(sampleOut => 1);
}
print "All done.\n" . $stats->Show();
