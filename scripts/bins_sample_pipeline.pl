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

This script detects the type of the sample based on the fastq files in the directory.  There can be only one
fastq file of each type. The following formats are supported.

=over 4

=item HMP

Paired-end reads with singletons. I<*.1.fastq> for the left reads, I<*.2.fastq> for the right reads, and
I<*.singleton.fastq> for the singletons.

=item MH

Paired-end reads. I<*_1.fastq> for the left reads and I<*_2.fastq> for the right reads.

=item SYNTH

Paired-end reads. C<synth1.fq> for the left reads and C<synth2.fq> for the right reads.

=item AG

Interleaved paired-end reads. C<*.fastq> for all of them.

=item PQ

Interleaved paired-end reads. C<*.fq> for all of them.

=item CT

Paired-end reads. C<*_R1_001.fastq> for the left reads and C<*_R2_001.fastq> for the right reads.

=back

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

=item reset

Delete all files except the assembly results (C<contigs.fasta> and C<output.contigs2reads.txt>) to force re-binning.
If C<all> is specified, even the assembly results will be deleted to force re-assembly, too.

=item gz

If specified, then the reads are stored in C<gz> files and must be unzipped first.

=item engine

Type of binning engine to use-- C<s> for standard or C<2> for alternate.

=item noIndex

If specified, the annotated genomes will not be indexed in PATRIC.

=item noBin

If specified, samples will be assembled but not binned.

=item nameSuffix

The suffix to use when naming the bins:  this is appended to the species name.  The default is C<clonal population>.

=item seedrole

The name of the universal role used for seeding the bins (generally only needed for status messages). The default is
C<PhenTrnaSyntAlph>.

=item seedProtFasta

A FASTA file containing examples of the universal role to use for seeding the bin assignment.

=item seedfasta

The name of the BLAST database for the seed protein in the various PATRIC genomes. The default is
C<PhenTrnaSyntAlph.fa> in the global data directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleID workDir',
        ["user|u=s", "user name for RAST access", { default => $ENV{RASTUSER} }],
        ["password|p=s", "password for RAST access", { default => $ENV{RASTPASS} }],
        ["force", "rebuild all files"],
        ["project=s", "source project type", { default => 'MH' }],
        ["reset:s", "delete all files except the assembly results to force re-binning"],
        ["gz", "unzip read files before processing"],
        ['engine=s', 'type of binning engine to use', { default => 's' }],
        ['noIndex', 'do not index annotated genomes in PATRIC'],
        ['noBin', 'do not bin the samples'],
        ['nameSuffix=s', 'suffix to assign to the bin name', { default => 'clonal population' }],
        ['seedProtFasta=s', 'name of a FASTA file containing examples of the seed protein to use for seeding the bins',
                            { default => "$FIG_Config::p3data/seedProt.fa" }],
        ['seedrole|R=s',   'ID of the universal role to seed the bins', { default => 'PhenTrnaSyntAlph' }],
        ['seedfasta|F=s',  'BLAST database (or FASTA file) of seed protein in all genomes', { default => "$FIG_Config::p3data/PhenTrnaSyntAlph.fa"}],
        );
# Close STDIN so we don't hang a parent process.
close(STDIN);
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
# Are we unzipping?
if ($opt->gz) {
    # Yes. Get the list of gz files.
    opendir(my $dh, $workDir) || die "Could not open work directory: $!";
    my @files = grep { $_ =~ /\.gz$/ } readdir $dh;
    my $total = scalar @files;
    my $count = 1;
    for my $file (@files) {
        print "Decompressing $file ($count of $total).\n";
        my $rc = system('gunzip', "$workDir/$file");
        die "Gunzip failed for $file: rc = $rc." if $rc;
        $count++;
    }
}
# Compute the FASTQ file names based on guesses about the project.
# We only need to do this if there are not yet assembled contigs.
if (! -s "$workDir/contigs.fasta") {
    SamplePipeline::PrepareAssembly($workDir, \%options);
}
# Store the engine.
$options{engine} = $opt->engine;
# Store the options.
$options{noIndex} = $opt->noindex // 0;
$options{noBin} = $opt->nobin // 0;
$options{nameSuffix} = $opt->namesuffix;
$options{seedProtFasta} = $opt->seedprotfasta;
$options{seedfasta} = $opt->seedfasta;
$options{seedrole} = $opt->seedrole;
# Check the file names. Save the file lengths.
my $resetOpt = $opt->reset;
# Are we resetting?
if (defined $resetOpt) {
    # We don't use this file any more, but we save it from deletions anyway.
    my $expectF = "$workDir/$sampleID" . "_abundance_table.tsv";
    # Yes. Get the list of files and delete the binning stuff.
    opendir(my $dh, $workDir) || die "Could not open work directory: $!";
    my @files = grep { -f "$workDir/$_" } readdir $dh;
    print "Deleting intermediate files in $workDir.\n";
    my ($count, $total) = (0,0);
    for my $file (@files) {
        my $fullName = "$workDir/$file";
        $total++;
        unless ($fullName eq $expectF || $fullName =~ /\.fastq$/ || $fullName =~ /\.fq$/ || $file eq 'site.tbl' ||
                $file eq 'run.log' || $file eq 'err.log' || $file eq 'exclude.tbl') {
            if ($resetOpt || ($file ne 'contigs.fasta' && $file ne 'output.contigs2reads.txt')) {
                unlink $fullName;
                $count++;
            }
        }
    }
    print "$count of $total files deleted.\n";
}
# Process the pipeline.
SamplePipeline::Process($workDir, %options);