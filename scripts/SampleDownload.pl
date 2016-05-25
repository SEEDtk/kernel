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
use File::Spec;

=head1 Download Human Microbiome Samples

    SampleDownload.pl [ options ] workDir

This script downloads Human Microbiome Project samples from the NCBI web site. It accepts as input
a tab-delimited file consisting of a site ID in the first column (e.g. C<palatine_tonsils>) and a
sample ID in the second column (e.g. C<SRS015061>) for each sample to be downloaded. The samples
are downloaded and unpacked. They must then be fed through the pipeline individually. (This will
presumably be done in parallel.)

=head2 Parameters

The positional parameter is the name of a work directory. Each sample is loaded into a subdirectory of this 
work directory. The command-line options are those found in L<ScriptUtils/ih_options> (for the standard input)
plus the following.

=over 4

=item force

If specified, samples will be downloaded even if their directories already exist.

=back 

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir', ScriptUtils::ih_options(),
                            ['force', 'force downloading over existing directories'],
        );
# Get the force option.
my $force = $opt->force;
# Verify the work directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir specified.";
}
# Get an absolute path for the working directory.
$workDir = File::Spec->rel2abs($workDir);
print "Absolute working directory path is $workDir.\n";
# Count the samples.
my $count = 0;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# These are our options for curl.
my @curlOpts = ('--remote-name', '--silent', '--show-error');
# Loop through the input.
while (! eof $ih) {
    # Switch to the working directory.
    chdir $workDir;
    # Get the next sample.
    my $line = <$ih>;
    chomp $line;
    my ($site, $sample) = split /\t/, $line;
    $count++;
    print "Processing $count: $sample from $site.\n";
    my $sampleDir = "$workDir/$sample";
    my $download = $force;
    if (! -d $sampleDir) {
        File::Copy::Recursive::pathmk($sampleDir);
        $download = 1;
    }
   # Now we have a directory into which we can put the files.
    if (! $download) {
        # We can skip this one.
        print "Directory skipped.\n";
    } else {
        # Downloading is required. Download and unpack the sample reads.
        print "Downloading samples.\n";
        my $rc = system('curl', @curlOpts, "http://downloads.hmpdacc.org/data/Illumina/$site/$sample.tar.bz2");
        die "Error code $rc downloading sample." if $rc;
        $rc = system('tar', '-xjvf', "$sample.tar.bz2");
        die "Error code $rc unpacking sample." if $rc;
        # Delete the tar file.
        unlink "$sample.tar.bz2";
        # Create the site file.
        print "Creating site file.\n";
        open(my $oh, '>', "$sample/site.tbl") || die "Could not open site file: $!";
        my $siteName = join(' ', map { ucfirst $_ } split /_/, $site);
        print $oh "HMP\t$site\t$siteName\n";
        close $oh;
        # Download the abundance profile.
        print "Downloading profile.\n";
        chdir $sampleDir;
        my $abundanceF = $sample . '_abundance_table.tsv.bz2';
        $rc = system('curl', @curlOpts, "http://downloads.hmpdacc.org/data/HMSCP/$sample/$abundanceF");
        die "Error code $rc downloading abundance." if $rc;
        $rc = system('tar', '-xjvf', $abundanceF);
        die "Error code $rc downloading abundance." if $rc;
        # Delete the tar file.
        unlink $abundanceF;
    }
}
print "All done.\n";