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

=item project

The name of the project. The default is C<HMP> for the Human Microbiome project. If C<MH> is specified, then
the samples will be presumed to be from the MetaHit project, which has a very different format.

=back 

=cut

use constant PROJECTS => { 'HMP' => 1, 'MH' => 1 };

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir', ScriptUtils::ih_options(),
                            ['force', 'force downloading over existing directories'],
                            ['project=s', 'ID of the source project', {default => 'HMP'}],
        );
# Get the force option.
my $force = $opt->force;
# Determine the project.
my $project = $opt->project;
if (! PROJECTS->{$project}) {
    die "Project $project not understood.";
} else {
    print "Project $project selected.\n";
}
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
        # We must create the directory in this case.
        File::Copy::Recursive::pathmk($sampleDir);
        $download = 1;
    } else {
        # Clear the current contents to avoid download/rename/unzip errors.
        File::Copy::Recursive::pathempty($sampleDir);
    }
   # Now we have a directory into which we can put the files.
    if (! $download) {
        # We can skip this one.
        print "Directory skipped.\n";
    } else {
        # Downloading is required. Download and unpack the sample reads.
        if ($project eq 'HMP') {
            print "Downloading samples.\n";
            my $rc = system('curl', @curlOpts, "http://downloads.hmpdacc.org/data/Illumina/$site/$sample.tar.bz2");
            die "Error code $rc downloading sample." if $rc;
            print "Unpacking samples.\n";
            $rc = system('tar', '-xjvf', "$sample.tar.bz2");
            die "Error code $rc unpacking sample." if $rc;
            # Delete the tar file.
            unlink "$sample.tar.bz2";
        } elsif ($project eq 'MH') {
            # For MetaHit, we need to compute the subdirectory from the project ID.
            my $subdir = substr($sample, 0, 6);
            # Now process the two samples. They are in separate files.
            for my $type (qw(1 2)) {
                print "Downloading sample $type.\n";
                my $fileName = $sample . "_$type.fastq";
                my $rc = system('curl', @curlOpts, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$subdir/$sample/$fileName.gz");
                die "Error code $rc downloading sample $type." if $rc;
                print "Unpacking sample $type.\n";
                # Note this deletes the gz file automatically.
                $rc = system('gunzip', "$fileName.gz");
                die "Error code $rc unpacking sample $type." if $rc;
                # Convert to HMP naming conventions.
                rename($fileName, "$sample.denovo_duplicates_marked.trimmed.$type.fastq") ||
                    die "Failed to renamed sample $type: $!";
            }
        } else {
            die "Project $project not implemented for sample download.";
        }
        # Create the site file.
        print "Creating site file.\n";
        open(my $oh, '>', "$sample/site.tbl") || die "Could not open site file: $!";
        my $siteName = join(' ', map { ucfirst $_ } split /_/, $site);
        print $oh "$project\t$site\t$siteName\n";
        close $oh;
        # Download the abundance profile (HMP only).
        if ($project eq 'HMP') {
            print "Downloading profile.\n";
            chdir $sampleDir;
            my $abundanceF = $sample . '_abundance_table.tsv.bz2';
            my $rc = system('curl', @curlOpts, "http://downloads.hmpdacc.org/data/HMSCP/$sample/$abundanceF");
            die "Error code $rc downloading abundance." if $rc;
            print "Unpacking abundance.\n";
            # Note this deletes the bz2 file automatically.
            $rc = system('bzip2', '--decompress', '--force', $abundanceF);
            die "Error code $rc unpacking abundance." if $rc;
        }
    }
}
print "All done.\n";