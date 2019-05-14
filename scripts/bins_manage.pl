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
use File::Copy::Recursive;


=head1 Process All Samples in a Binning Directory

    bins_manage.pl [ options ] binDir

This is a long-running job that processes the binning pipeline in a binning directory.  Every 5 minutes it will wake
up and check for progress, then start or restart jobs in order to keep the binning pipeline going.  When all the
samples are complete it will stop.  It can also be stopped using the C<--stop> option.

=head2 Parameters

The positional parameter is the name of the binning directory.  All the samples must be in subdirectories.

The command-line options are the following

=over 4

=item maxJobs

The maximum number of jobs to run at one time.  The default is C<20>.

=item maxAsm

The maximum number of assemblies to run at one time.  The default is C<4>.  Note that an assembly counts as a job,
so if this parameter is 4 and the maximum number of jobs is 20, there can be 16 non-assembly jobs when all 4
assemblies are running.

=item noIndex

If specified, the bins will not be indexed in PATRIC.

=item sleep

The number of minutes to sleep between status checks.  The default is C<5>.

=item stop

If specified, a C<STOP> file will be written to the directory and this program will exit.  This will cause any
other instance of this script to stop when it next wakes up.

=back

=cut

# Files to keep during a reset.
use constant KEEPERS => { 'site.tbl' => 1, 'run.log' => 1, 'err.log' => 1, 'exclude.tbl' => 1, 'contigs.fasta' => 1,
        'output.contigs2reads.txt' => 1 };

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir', Shrub::script_options(),
        ['maxJobs=i', 'maximum number of jobs to have running at one time', { default => 20 }],
        ['maxAsm=i', 'maximum number of assemblies to have running at one time'],
        ['noIndex', 'if specified, the bins produced will not be indexed in PATRIC'],
        ['sleep=i', 'number of minutes to sleep between awakening cycles'],
        ['stop', 'attempt to stop another management process']
        );
# Create a statistics object.
my $stats = Stats->new();
# Get the input directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir missing or invalid.";
}
# Process the stop option first.
if ($opt->stop) {
    open (my $oh, '>', "$binDir/STOP") || die "Could not open stop file: $!";
    print $oh "STOP\n";
    print "STOP file created.\n";
} else {
    # Here we want to go into the loop.  Get the options.
    my $maxJobs = $opt->maxjobs;
    my $maxAsm = $opt->maxasm;
    my $noIndex = ($opt->noindex ? '--noIndex' : '');
    my $sleep = $opt->sleep * 60;
    # Loop until we find the stop file.
    while (! -f "$binDir/STOP") {
        $stats->Add(cycles => 1);
        # This will count the incomplete directories.
        my $incomplete = 0;
        # Get a hash of the running jobs.
        my %running;
        my $jobsLeft = $maxJobs;
        my $asmLeft = $maxAsm;
        my @jobs = `ps -AF`;
        for my $job (@jobs) {
            if ($job =~ /bins_sample_pipeline\s+(?:--\S+\s+)*(\w+)/) {
                my $sample = $1;
                $running{$sample} = 1;
                $jobsLeft--;
                $incomplete++;
                # Is an assembly in progress?
                if (-d "$binDir/$sample/Assembly" && ! -s "$binDir/$sample/contigs.fasta") {
                    $asmLeft--;
                }
            }
        }
        # Only proceed if we have room to start a job.
        if ($jobsLeft > 0) {
            # Get the subdirectories.
            opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
            my @samples = sort grep { substr($_,0,1) ne '.' && -d "$binDir/$_" } readdir $dh;
            closedir $dh;
            # We now run through the directories.  Directories that are still running are skipped.  If a directory is
            # in the downloaded state, we queue it for assembly. If it is evaluating, we reset it and queue it for
            # resume.  If it is assembled and not done, we simply queue it for resume.
            my (@assemble, @resume);
            for my $sample (@samples) {
                my $subDir = "$binDir/$sample";
                # Only process the directory if it is not running.
                if (! $running{$sample}) {
                    if (-s "$subDir/contigs.fasta" && -d "$subDir/Assembly") {
                        # We have a contig file and there is assembly data.  Clean it up.
                        ClearAssembly($subDir);
                        opendir(my $dh, $subDir) || die "Could not open $subDir: $!";
                        my @fastqs = grep { $_ =~ /\.(?:fastq|fq)$/ } readdir $dh;
                        for my $fastq (@fastqs) {
                            unlink "$subDir/$fastq";
                        }
                        $stats->Add(dirsCleaned => 1);
                        print "$sample reads cleaned.\n";
                    }
                    # At this point, we only care about the directory if it is incomplete.
                    if (! -s "$subDir/Eval/index.tbl") {
                        $incomplete++;
                        if (-f "$subDir/bins.rast.json") {
                            # Here RAST is complete, but we failed during evaluation.  This means we have to rebin.
                            print "Must rebin $sample.\n";
                            # Delete the assembly work directory.
                            ClearAssembly($subDir);
                            # Get the list of files and delete the binning stuff.
                            opendir(my $dh, $subDir) || die "Could not open work directory: $!";
                            my @files = grep { -f "$subDir/$_" } readdir $dh;
                            for my $file (@files) {
                                my $fullName = "$subDir/$file";
                                unless ($fullName =~ /_abundance_table.tsv$/ || $fullName =~ /\.fastq$/ || $fullName =~ /\.fq/ ||
                                        KEEPERS->{$file}) {
                                    unlink $fullName;
                                }
                            }
                            # Queue for resume.
                            push @resume, $sample;
                        } elsif (! -s "$subDir/contigs.fasta" && -s "$subDir/site.tbl") {
                            # This directory is unassembled.
                            if (-d "$subDir/Assembly") {
                                # The assembly crashed, so prep for restart.
                                print "Must backout $sample.\n";
                                ClearAssembly($subDir);
                                $stats->Add(assemblyBackout => 1);
                            }
                            # Queue for assembly.
                            push @assemble, $sample;
                        } else {
                            # Queue for resume.
                            push @resume, $sample;
                        }
                    }
                }
            }
            while ($jobsLeft && $asmLeft && @assemble) {
                my $sample = shift @assemble;
                StartJob($binDir, $sample, $noIndex, 'Started');
                $jobsLeft--;
                $asmLeft--;
            }
            while ($jobsLeft && @resume) {
                my $sample = shift @resume;
                StartJob($binDir, $sample, $noIndex, 'Restarted');
                $jobsLeft--;
            }
        }
        # Wait for the next wakeup.
        sleep $sleep;
    }
}
print "All done.\n" . $stats->Show();

sub StartJob {
    my ($binDir, $dir, $noIndex, $start) = @_;
    my $subDir = "$binDir/$dir";
    my $cmd = "bins_sample_pipeline $noIndex $dir $subDir >$subDir/run.log 2>$subDir/err.log";
    my $rc = system("nohup $cmd &");
    print "$start job for $dir.\n";
    $stats->Add("jobs$start" => 1);
}


sub ClearAssembly {
    my ($subDir) = @_;
    File::Copy::Recursive::pathempty("$subDir/Assembly") || die "Could not empty $subDir/Assembly: $!";
    rmdir "$subDir/Assembly";
}