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
use JSON::XS;
use Cwd 'abs_path';
use GPUtils;

=head1 Produce a Quality Report on Genome Packages

    package_report.pl [ options ] dir package

This script examines the information in a genome package to produce a one-line quality report. The report is then stored in the
C<quality.tbl> file inside the package directory. If the file exists, it will simply be re-used (unless the C<--force> option is
specified). If no specific package ID is specified, the entire package directory is examined and a multi-line report is
produced.

The columns in the report are as follows.

=over 4

=item 1

Sample name (if any)

=item 2

Genome ID

=item 3

Genome name

=item 4

Number of contigs

=item 5

Number of DNA base pairs

=item 6

C<1> if the genome is mostly complete, else C<0>

=item 7

The contig N50 score.

=item 8

The contig N70 score.

=item 9

The contig N90 score.

=item 10

ID of closest reference genome

=item 11

Name of closest reference genome

=item 12

SciKit coarse evaluation score (percent)

=item 13

SciKit fine evaluation score (percent)

=item 14

CheckG evaluation completeness (percent)

=item 15

CheckG evaluation contamination (percent)

=item 16

CheckG taxonomy classification

=item 17

C<1> for a good genome, C<0> for a bad one.

=back

The genome ID is the same as the subdirectory (package) name. The other values are taken from the following files.

=over 4

=item data.tbl

A tab-delimited file of key/value pairs. The following keys are important.

=over 8

=item *

Sample Name

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

=item EvalByCheckG/evaluate.log

CheckG results file. The second line of this file contains tab-delimited fields, containing (0) the completeness
score, (1) the contamination score, and (3) the taxonomic grouping used. Note that the third field (2) is not used.

=item EvalBySciKit/evaluate.log

SciKit evaluator results file, tab-delimited. At the end of the file is a series of lines containing scores. The line beginning with
C<Coarse_Consistency> contains the coarse score and the line beginning with C<Fine_Consistency> contains the fine score.

=item bin.gto

A L<GenomeTypeObject> in json format, used to compute completeness and the N-metrics.

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

=item samples

If specified, a comma-delimited list of sample names. Only packages from the specified sample will be output.

=item json

If specified, emit a JSON-formatted report instead of text.

=item log

If specified, progress will be written to the standard error output.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dir package',
        ["force", 'force regeneration of quality data'],
        ["quiet", 'suppress standard output'],
        ["missing", 'list missing data to error output'],
        ["samples=s", 'only process packages from the specified samples'],
        ["json", 'write the report in json format'],
        ["log|verbose|v", "show progress on STDERR"],
        );
# Get the log option.
my $log = $opt->log;
# Get the directory and the package.
my $json_packages = [];
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
        print STDERR "Reading directory $dir.\n" if $log;
        # Here we want multiple packages.
        opendir(my $dh, $dir) || die "Could not open package directory: $!";
        # Check for filtering.
        if (! $opt->samples) {
            # No filtering. Get all the packages.
            @packages = sort grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
        } else {
            # Loop through the packages, seeking ones from the specified sample.
            print STDERR "Filtering for samples.\n" if $log;
            my %samples = map { $_ => 1 } split /,/, $opt->samples;
            for my $package (sort grep { $_ =~ /^\d+\.\d+$/ } readdir $dh) {
                if (open(my $ih, '<', "$dir/$package/data.tbl")) {
                    while (! eof $ih) {
                        my $line = <$ih>;
                        if ($line =~ /^Sample Name\t(.+)/) {
                            if ($samples{$1}) {
                                push @packages, $package;
                            }
                        }
                    }
                }
            }
        }
    }
    # Loop through the packages, producing output.
    my $total = scalar @packages;
    my $count = 0;
    print STDERR "$total packages to process.\n" if $log;
    for my $package (@packages) {
        my $line = produce_report($dir, $package, $force, $report);
        if ($opt->json)	{
            push(@$json_packages, $line);
        } else {
            if (! $opt->quiet) {
                print $line;
            }
        }
        $count++;
        print STDERR "$count of $total packages processed.\n" if $log && ($count % 100 == 0);
    }
}
if ($opt->json) {
    my $path = $dir;
    $path = abs_path($path) unless $path =~ m,^/,;
    my $json_report = {
        package_directory => abs_path($path),
        packages => $json_packages,
    };
    print JSON::XS->new->pretty(1)->canonical(1)->encode($json_report);
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
            # Compute the GTO metrics.
            my $gto = GenomeTypeObject->create_from_file("$pDir/bin.gto");
            my $metricsH = $gto->metrics();
            # This will list the missing scores.
            my @missing;
            # Get the checkm scores.
            my ($checkGscore, $checkGcontam, $checkGtaxon) = ('', '', '');
            if (-d "$pDir/EvalByCheckG") {
                my $found;
                if (open(my $fh, '<', "$pDir/EvalByCheckG/evaluate.log")) {
                    my $line = <$fh>;
                    $line = <$fh>; chomp $line;
                    ($checkGscore, $checkGcontam, undef, $checkGtaxon) = split /\t/, $line;
                    $found = 1;
                    close $fh;
                }
                if (! $found) {
                    push @missing, 'CheckG';
                }
            }
            # Get the scikit score.
            my ($scikitCScore, $scikitFScore) = ('', '');
            if (-d "$pDir/EvalBySciKit") {
                my $found;
                if (open(my $fh, '<', "$pDir/EvalBySciKit/evaluate.log")) {
                    while (! eof $fh) {
                        my $line = <$fh>;
                        if ($line =~ /^Coarse_Consistency=\s+(.+)\%/) {
                            $scikitCScore = $1;
                            $found = 1;
                        } elsif ($line =~ /^Fine_Consistency=\s+(.+)%/) {
                            $scikitFScore = $1;
                        }
                    }
                    close $fh;
                }
                if (! $found) {
                    push @missing, 'SciKit';
                }
            }
            # Is this a good genomes?
            my $seedFlag = GPUtils::good_seed($gto);
            my $goodFlag = (($seedFlag && $scikitFScore && $scikitFScore >= 85 && $checkGscore && $checkGscore >= 80 && $checkGcontam <= 15) ? 1 : 0);
            # Assemble the output line.
            my $refGenome = $dataVals{'Ref Genome'} // $dataVals{'Source Package'} // $dataVals{'Source Database'} // '';
            my $refName = $dataVals{'Ref Name'} // 'derived';
            my $sampleName = $dataVals{'Sample Name'} // 'Derived';
            $retVal = join("\t", $sampleName, $package, $dataVals{'Genome Name'}, $dataVals{'Contigs'},
                    $dataVals{'Base pairs'}, $metricsH->{complete}, $metricsH->{N50}, $metricsH->{N70},
                    $metricsH->{N90}, $refGenome, $refName, $scikitCScore, $scikitFScore, $checkGscore, $checkGcontam,
                    $checkGtaxon, $goodFlag) . "\n";
            open(my $oh, '>', "$pDir/quality.tbl") || die "Could not write to quality file for $package: $!";
            print $oh $retVal;
            # Check for the missing report.
            if (scalar(@missing) && $report) {
                print STDERR join("\t", $package, @missing) . "\n";
            }
        }
    }
    #
    # If we are using json, split and emit a structure. Otherwise return the line.
    #

    if ($opt->json)
    {
        my $l = $retVal;
        chomp $l;
        my($sample_name, $genome_id, $genome_name, $contigs, $dna, $complete,
           $n50, $n70, $n90, $ref_id, $ref_name,
           $scikit_coarse, $scikit_fine, $checkg_completeness, $checkg_contamination, $checkg_tax) = split(/\t/, $l);
        $retVal = {
            sample_name => $sample_name,
            genome_id => $genome_id,
            genome_name => $genome_name,
            contigs => 0 + $contigs,
            dna_bp => 0 + $dna,
            complete => $complete,
            n50 => 0 + $n50,
            n70 => 0 + $n70,
            n90 => 0 + $n90,
            scikit_coarse => 0.0 + $scikit_coarse,
            scikit_fine => 0.0 + $scikit_fine,
            checkg_completeness => 0.0 + $checkg_completeness,
            checkg_contamination => 0.0 + $checkg_contamination,
            checkg_taxonomy => $checkg_tax,
        };
    }

    return $retVal;
}
