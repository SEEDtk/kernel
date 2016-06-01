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
use SamplePipeline;

=head1 Process a Metagenome Sample

    bins_sample_pipeline.pl [ options ] sampleID workDir

This script runs a single sample through the binning pipeline. This includes assembling the reads, forming the bins,
applying RAST to each bin, and checking expectations. The bulk of the work is performed by L<SamplePipeline>.

This script assumes the sample comes from the Human Microbiome Project and conforms to its naming conventions. In the
future, command-line options will be added to support other formats.

=head2 Parameters

The positional parameters are the ID of the sample and the name of the working directory.

The following command-line options are supported.

=over 4

=item user

User name for RAST access. If omitted, the default is taken from the RASTUSER environment variable.

=item password

Password for RAST access. If omitted, the default is taken from the RASTPASS environment variable.

=item force

Force rebuilding of all files, even ones that already exist. Otherwise, if files exist in the directory, it will
be presumed a previous run failed in progress and it will be resumed.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleID workDir',
        ["user|u=s", "user name for RAST access", { default => $ENV{RASTUSER} }],
        ["password|p=s", "password for RAST access", { default => $ENV{RASTPASS} }],
        ["force", "rebuild all files"],
        );
# Get the sample name and work directory.
my ($sampleID, $workDir) = @ARGV;
if (! $sampleID) {
    die "No sample ID specified.";
} elsif (! $workDir) {
    die "No work directory specified.";
} elsif (! -d $workDir) {
    die "Invalid work directory $workDir.";
}
# This will contain our options for the pipeline.
my %options = (force => $opt->force, user => $opt->user, password => $opt->password);
# Compute the file names for the sample.
my ($f1q, $f2q, $fsq) = map { "$workDir/$sampleID.denovo_duplicates_marked.trimmed.$_.fastq" } qw(1 2 singleton);
my $expectF = "$workDir/$sampleID" . "_abundance_table.tsv";
# Check the files. Save the ones that exist as options.
if (-f $f1q) {
    $options{f1} = $f1q;
}
if (-f $f2q) {
    $options{f2} = $f2q;
}
if (-f $fsq) {
    $options{fs} = $fsq;
}
if (-f $expectF) {
    $options{expect} = $expectF;
}
# Process the pipeline.
SamplePipeline::Process($workDir, %options);