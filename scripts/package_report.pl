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

=head1 Produce a Quality Report on Genome Packages

    package_report.pl [ options ] dir package

This script examines the information in a genome package to produce a one-line quality report. The report is then stored in the
C<quality.tbl> file inside the package directory. If the file exists, it will simply be re-used (unless the C<--force> option is
specified). If no specific package ID is specified, the entire package directory is examined and a multi-line report is
produced.

The columns in the report are as follows.

=over 4

=item 1

Genome ID

=item 2

Genome name

=item 3

Number of contigs

=item 4

Number of DNA base pairs

=item 5

ID of closest reference genome

=item 6

Name of closest reference genome

=item 7

SciKit evaluation score (percent)

=item 8

Tensor Flow evaluation score (percent)

=item 9

CheckM evaluation completeness (percent)

=item 10

CheckM evaluation contamination (percent)

=item 11

CheckM taxonomy classification

=back

The genome ID is the same as the subdirectory (package) name. The other values are taken from the following files.

=over 4

=item data.tbl

A tab-delimited file of key/value pairs. The following keys are important.

=over 8

=item *

Genome Name

=item *

Contigs

=item *

Base pairs

=item *

Ref Genome

=item *

Ref Name

=back

=item EvalByCheckm/evaluate.log

CheckM results file. The fourth line of this file contains space-delimited fields, the first of which is the word C<bin>. The second field is the
CheckM taxonomy classification. The thirteenth is the checkM completeness and the fourteenth is the checkM contamination score.

=item EvalBySciKit/evaluate.log

SciKit evaluator results file, tab-delimited. At the end of the file is a series of lines containing scores. The line beginning with C<Pct_roles=>
C<Consistency=> contains the score.

=item EvalByTF/evaluate.log

Tensor Flow evaluator results file, containing a single number--the score.

=back

=head2 Parameters

The positional parameters are the name of the directory containing the genome packages and the optional package ID. If no package ID is
specified, all packages are processed.

The command-line options are as follows.

=over 4

=item force

If no single package is specified, then rather than generate quality data for each package, the data from existing C<quality.tbl> files is
used. If this option is specified, quality data is regenerated for everything.

=item quiet

If specified, no standard output is generated.

=item missing

If specified, a list of the directories with missing information will be sent to the standard error output.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir package',
        ["force", 'force regeneration of quality data'],
        ["quiet", 'suppress standard output'],
        ["missing", 'list missing data to error output']
        );
# Get the directory and the package.
my ($dir, $package) = @ARGV;
if (! $dir) {
    die "No packages directory specified.";
} elsif (! -d $dir) {
    die "Invalid package directory $dir.";
} else {
    # We will put the report lines in here.
    my @lines;
    # This is the force flag. We default to the option value.
    my $force = $opt->force;
    # This is the missing-item report flag.
    my $report = $opt->missing;
    # This will contain the list of packages to process.
    my @packages;
    if ($package) {
        # Here we have one package. Turn on forcing.
        push @packages, $package;
        $force = 1;
    } else {
        # Here we want all the packages.
        opendir(my $dh, $dir) || die "Could not open package directory: $!";
        @packages = grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
    }
    # Loop through the packages, producing output.
    for my $package (@packages) {
        my $line = produce_report($dir, $package, $force, $report);
        if (! $opt->quiet) {
            print $line;
        }
    }
}

# Produce the report for one package.
sub produce_report {
    my ($dir, $package, $force, $report) = @_;
    # This will be the return line.
    my $retVal;
    # Create the package directory name.
    my $pDir = "$dir/$package";
    if (! $force && open(my $ih, '<', "$pDir/quality.tbl")) {
        # Here we have a precomputed result.
        $retVal = <$ih>;
    } else {
        my $dataFile = "$pDir/data.tbl";
        if (! -s $dataFile) {
            die "Invalid package $package: missing data.tbl file.";
        } else {
            # Read the data.tbl file.
            open(my $ih, '<', $dataFile) || die "Could not open data file for $package: $!";
            my %dataVals;
            while (! eof $ih) {
                my $line = <$ih>;
                chomp $line;
                my ($key, $value) = split /\t/, $line;
                $dataVals{$key} = $value;
            }
            # This will list the missing scores.
            my @missing;
            # Get the checkm scores.
            my ($checkMscore, $checkMcontam, $checkMtaxon) = ('', '', '');
            if (-d "$pDir/EvalByCheckm") {
                my $found;
                if (open(my $fh, '<', "$pDir/EvalByCheckm/evaluate.log")) {
                    while (! eof $fh) {
                        my $line = <$fh>;
                        if ($line =~ /^\s+bin\s+/) {
                            my @cols = split /\s+/, $line;
                            $checkMtaxon = $cols[2];
                            $checkMscore = $cols[13];
                            $checkMcontam = $cols[14];
                            $found = 1;
                        }
                    }
                    close $fh;
                }
                if (! $found) {
                    push @missing, 'Checkm';
                }
            }
            # Get the scikit score.
            my $scikitScore = '';
            if (-d "$pDir/EvalBySciKit") {
                my $found;
                if (open(my $fh, '<', "$pDir/EvalBySciKit/evaluate.log")) {
                    while (! eof $fh) {
                        my $line = <$fh>;
                        if ($line =~ /^Consistency=\s+(.+)\%/) {
                            $scikitScore = $1;
                            $found = 1;
                        }
                    }
                    close $fh;
                }
                if (! $found) {
                    push @missing, 'SciKit';
                }
            }
            # Get the tensor flow score.
            my $tfScore = '';
            if (-d "$pDir/EvalByTF") {
                my $found;
                if (open(my $fh, '<', "$pDir/EvalByTF/evaluate.log")) {
                    my $line = <$fh>;
                    chomp $line;
                    $tfScore = $line;
                    close $fh;
                    $found = 1;
                }
                if (! $found) {
                    push @missing, 'TF';
                }
            }
            # Assemble the output line.
            my $refGenome = $dataVals{'Ref Genome'} // $dataVals{'Source Package'} // $dataVals{'Source Database'} // '';
            my $refName = $dataVals{'Ref Name'} // 'derived';
            $retVal = join("\t", $package, $dataVals{'Genome Name'}, $dataVals{'Contigs'}, $dataVals{'Base pairs'},
                    $refGenome, $refName, $scikitScore, $tfScore, $checkMscore, $checkMcontam, $checkMtaxon) .
                    "\n";
            open(my $oh, '>', "$pDir/quality.tbl") || die "Could not write to quality file for $package: $!";
            print $oh $retVal;
            # Check for the missing report.
            if (scalar(@missing) && $report) {
                print STDERR join("\t", $package, @missing) . "\n";
            }
        }
    }
    # Return the output line.
    return $retVal;
}
