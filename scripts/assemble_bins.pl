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
use SeedAware;
use File::Copy::Recursive;

=head1 Assemble Bins from Reads

    assemble_bins.pl [ options ] binDir

This script will assemble the bins output from L<approx_bin_DNA.pl>. The bins are in interlaced FASTQ files in
the specified bin directory. This script will loop through the files, assembling and creating contig files with
the same name but the suffix C<.fasta> instead of C<.fastq>.

=head2 Parameters

The positional parameter is the name of the directory containing the bin FASTQ files. All files in this directory
with the suffix C<.fastq> will be processed. A working sub-directory with the name C<Assemble> will be created
or used as the SPAdes working directory. The output contigs will be produced in a file named C<contigs.fasta> in
this directory and copied to the output directory before the working sub-directory is cleared.

The command-line options are as follows.

=over 4

=item resume

Only process a FASTQ file if no matching FASTA file exists.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ['resume', 'assemble unprocessed input files only'],
        );
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "Invalid or missing input directory $binDir.";
}
# Create the working directory.
my $workDir = "$binDir/Assemble";
if (! -d $workDir) {
    print "Creating working directory $workDir.\n";
    File::Copy::Recursive::pathmk($workDir);
}
# Get the input files.
opendir(my $dh, $binDir) || die "Could not open input directory $binDir: $!";
my @inputs = grep { $_ =~ /\.fastq$/ } readdir $dh;
closedir $dh;
print scalar(@inputs) . " bin files found in $binDir.\n";
# Find the SPAdes assembler.
my $cmdPath = SeedAware::executable_for('spades.py');
print "Assembler found at $cmdPath.\n";
# Loop through the input files.
for my $inFile (@inputs) {
    # Compute the output file name.
    my $outFile = $inFile;
    substr($outFile, -1, 1) = 'a';
    # Create the parameters.
    my @parms = ('-o', $workDir, '--12', "$binDir/$inFile");
    # Call the assembler.
    print "Invoking the SPAdes assembler.\n";
    my $rc = system($cmdPath, @parms);
    die "Error exit $rc from SPAdes." if $rc;
    # Copy the output file.
    print "Cleaning up from assembly.\n";
    File::Copy::Recursive::fmove("$workDir/contigs.fasta", "$binDir/$outFile") || die "Could not copy assembly to $outFile: $!";
    # Clear the working directory.
    File::Copy::Recursive::pathempty($workDir);
}
print "All done.\n";