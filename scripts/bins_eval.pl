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
use Bin;
use File::Copy::Recursive;

=head1 Evaluate Bins After Downloading

    bins_eval.pl [ options ] binDir jobID group

This script runs the new evaluator on downloaded genome files and then runs them all through the summary-page application.
Between L</bins_fasta.pl> and this script, the output genomes from the annotation jobs need to have been downloaded from
BV-BRC, which is a long, tedious process.

=head2 Parameters

The positional parameters are the name of the directory containing the binning results (including the C<.genome> files), the
job ID to use, and the optional group name.

The command-line options are the following

=over 4

=item help

Display command-line usage.

=item eval

Name of the evaluation directory containing the evaluator files.  The default is C<Eval> in the PATRIC shared data directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir jobID groupName',
        ['--eval=s', 'evaluator directory', { default => "$FIG_Config::p3data/Eval" } ],
        );
# Get the parameters.
my ($binDir, $jobID, $groupName) = @ARGV;
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir is not found or is invalid.";
}
# Get the evaluation directory.
my $eval = $opt->eval;
if (! -d $eval) {
    die "Evaluation directory is not found or invalid.";
}
# Insure we have a work directory.
my $workDir = "$binDir/QTemp";
if (! -d $workDir) {
    File::Copy::Recursive::pathmk($workDir) || die "Could not make work directory: $!";
}
if (-s "$workDir/eval.log") {
    unlink "$workDir/eval.log";
}
my $outDir = "$binDir/QEval";
if (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) || die "Could not make output directory: $!";
}
# Read in the bins and build a hash of the reference genome IDs.
my $refHash = readBins($binDir);
# Get the list of genome files.
opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
my @binFiles = grep { $_ =~ /^Bin\.\d+\.\d+\.genome$/ } readdir $dh;
closedir $dh; undef $dh;
# Start the control file for the summary.
open(my $oh, '>', "$workDir/control.tbl") || die "Could not open output control file: $!";
print $oh "genomeFile\turl\n";
# Loop through the files, processing each one.
for my $binFile (@binFiles) {
    my $fileName = "$binDir/$binFile";
    if ($binFile =~ /^Bin\.(\d+)\.(\d+)/) {
        # Get the useful parts of the bin name.
        my $binNum = $1;
        my $taxID = $2;
        # Compute the ref ID and URL.
        my $refID = $refHash->{$taxID};
        my $url = "bin.$binNum.$taxID.html";
        my $gtoFile = "$outDir/bin.$binNum.$taxID.gto";
        print STDOUT "Processing bin $binNum with tax ID $taxID and reference genome $refID.\n";
        # Call the evaluator.
        my $rc = system("p3x-eval-gto", "--workDir", $workDir, "--ref", $refID, "--deep", "--evalDir", $eval,
                "--bins", "$binDir/bins.json", "--improve", $fileName, $gtoFile, "$outDir/$url");
        if ($rc) {
            die "Error in evaluation.  rc = $rc.";
        }
        # Write the control record.
        print $oh "$gtoFile\t$url\n";
    }
}
# Close the control file.
close $oh; undef $oh;
# Call the summary report.
my @parms = ("--job", $jobID);
if ($groupName) {
    push @parms, "--group", $groupName;
}
push @parms, "-i", "$workDir/control.tbl", "-o", "$outDir/index.html";
print STDERR "Running evaluator: binReport " . join(" ", @parms) . "\n";
my $rc = system("dl4j.eval", "binReport", @parms);
print STDERR "rc = $rc\n";

# Read in the binning directory and return a hash reference mapping taxon IDs to reference genomes.
sub readBins {
    my ($binDir) = @_;
    my %retVal;
    print STDOUT "Reading bins.json from $binDir.\n";
    my $bins = Bin::ReadBins("$binDir/bins.json");
    for my $bin (@$bins) {
        my $taxID = $bin->taxonID;
        my @refs = $bin->refGenomes;
        $retVal{$taxID} = $refs[0];
    }
    return \%retVal;
}
