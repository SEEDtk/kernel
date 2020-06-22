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
use Math::Round;
use SamplePipeline;

=head1 Process All Samples in a Binning Directory

    bins_manage.pl [ options ] binDir

This is a long-running job that processes the binning pipeline in a binning directory.  Every 2 minutes it will wake
up and check for progress, then start or restart jobs in order to keep the binning pipeline going.  When all the
samples are complete it will stop.  It can also be stopped using the C<--stop> option in L<bins_status.pl>.

=head2 Parameters

The positional parameter is the name of the binning directory.  All the samples must be in subdirectories.

The command-line options are the following

=over 4

=item maxJobs

The maximum number of jobs to run at one time.  The default is C<20>.

=item maxAsm

The maximum number of assemblies to run at one time.  The default is C<2>.  Note that an assembly counts as a job,
so if this parameter is 2 and the maximum number of jobs is 20, there can be 18 non-assembly jobs when both
assemblies are running.

=item noIndex

If specified, the bins will not be indexed in PATRIC.

=item sleep

The number of minutes to sleep between status checks.  The default is C<2>.

=item altseed

Name of the alternate seed protein to use.  If specified, there must be a DNA fasta file named I<protName>C<.fa>
and a corresponding protein fasta file named I<protName>C<.faa> in the SEEDtk global directory.  These files
will be used to seed the binning process and find reference genomes.

=back

=cut

# Files to keep during a reset.
use constant KEEPERS => { 'site.tbl' => 1, 'run.log' => 1, 'err.log' => 1, 'exclude.tbl' => 1, 'contigs.fasta' => 1,
        'output.contigs2reads.txt' => 1 };

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ['maxJobs=i', 'maximum number of jobs to have running at one time', { default => 20 }],
        ['maxAsm=i', 'maximum number of assemblies to have running at one time', { default => 1 }],
        ['noIndex', 'if specified, the bins produced will not be indexed in PATRIC'],
        ['sleep=i', 'number of minutes to sleep between awakening cycles', { default => 2 }],
        ['altseed=s', 'ID of an alternate seed protein to use']
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
# Here we want to go into the loop.  Get the options.
my $maxJobs = $opt->maxjobs;
my $maxAsm = $opt->maxasm;
my $extras = ($opt->noindex ? '--noIndex' : '');
if ($opt->altseed) {
    my $altSeed = $opt->altseed;
    $extras .= " --seedrole=$altSeed --seedProtFasta=$FIG_Config::global/$altSeed.faa --seedfasta=$FIG_Config::global/$altSeed.fa";
}
my $sleep = $opt->sleep * 60;
print "Starting main loop for $maxJobs jobs with up to $maxAsm assemblies.\n";
# Loop until we find the stop file.
while (! -f "$binDir/STOP") {
    $stats->Add(cycles => 1);
    # This will count the incomplete directories.
    my $incomplete = 0;
    # Get a hash of the running jobs.
    my %running;
    my $jobsLeft = $maxJobs;
    my $asmLeft = $maxAsm;
    my @jobs = `ps -Af`;
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
        # This will track the number of completed samples in each category.
        my %cats;
        # This will contain the category of each incomplete sample.
        my %sampCats;
        # We now run through the directories.  Directories that are still running are skipped.  If a directory is
        # in the downloaded state, we queue it for assembly. If it is evaluating, we reset it and queue it for
        # resume.  If it is binned and not done, we queue it for resume.  If it is assembled and not binned, we
        # queue it for startup.
        my (@assemble, @resume, @startup);
        for my $sample (@samples) {
            my $subDir = "$binDir/$sample";
            # Only process the directory if it is not running either here or on an assembly machine.
            if (! $running{$sample} && ! -f "$subDir/ASSEMBLE") {
                if (-s "$subDir/contigs.fasta" && -d "$subDir/Assembly") {
                    # We have a contig file and there is assembly data.  Clean it up.
                    SamplePipeline::ClearAssembly($subDir, $sample);
                    $stats->Add(dirsCleaned => 1);
                }
                # Get the sample category.
                my $site = "Unspecified";
                if (open(my $sh, '<', "$subDir/site.tbl")) {
                    my $line = <$sh>;
                    if ($line =~ /\t(\S+)\t/) {
                        $site = $1;
                    }
                }
                $cats{$site} //= 0;
                $sampCats{$sample} = $site;
                # We count the complete samples by category.
                if (-s "$subDir/Eval/index.tbl") {
                    # Here we are complete.
                    $cats{$site}++;
                    if (-f "$subDir/START") {
                        # Here it is the first time we've seen it completed, so write a message.
                        my $duration = '';
                        my @stats = stat("$subDir/START");
                        if (@stats) {
                            $duration = " in  " . Math::Round::nearest(0.1, (time - $stats[9])/3600) . " hours";
                        } else {
                            $duration = '';
                        }
                        if (! -s "$subDir/bin1.gto") {
                            $duration .= ' with no bins found';
                        }
                        print "$sample ($site) completed$duration.\n";
                        unlink "$subDir/START";
                    }
                } else {
                    # Here we are incomplete.  We need to check for a need to start this sample.
                    $incomplete++;
                    if (! -s "$subDir/site.tbl") {
                        # Here we are still downloading.
                    } elsif (-f "$subDir/bins.rast.json") {
                        # Here RAST is complete, but we failed during evaluation.  bins_status must fix this.
                    } elsif (! -s "$subDir/contigs.fasta") {
                        # This directory is unassembled.
                        if (-d "$subDir/Assembly") {
                            # The assembly crashed or it is running on an assembly machine.
                        } else {
                            # Queue for assembly.
                            push @assemble, $sample;
                        }
                    } elsif (-s "$subDir/bins.json") {
                        # Queue for resume. We can start annotating, so this is our highest priority.
                        push @resume, $sample;
                    } else {
                        # Queue for startup.
                        push @startup, $sample;
                    }
                }
            }
        }
        if ($jobsLeft) {
            # Here we need to start some jobs.  Do we have room for assemblies?
            if ($asmLeft) {
                # Yes.  Sort the samples from the rarest categories last.  We want to start those first.
                @assemble = sort { $cats{$sampCats{$b}} <=> $cats{$sampCats{$a}} } @assemble;
                # Start the jobs.
                while ($asmLeft > 0 && @assemble) {
                    my $sample = pop @assemble;
                    StartJob($binDir, $sample, $extras);
                    $asmLeft--;
                    $jobsLeft--;
                }
            }
            # Resume anything we have room for.
            push @resume, @startup;
            while ($jobsLeft > 0 && @resume) {
                my $sample = shift @resume;
                StartJob($binDir, $sample, $extras);
                $jobsLeft--;
            }
        }
        # Do we have any incomplete samples?
        if (! $incomplete) {
            # No, stop the loop.
            print "No samples left to process.  Stopping main loop.\n";
            StopFile();
        }
    }
    # Wait for the next wakeup.
    sleep $sleep;
}
print "Deleting stop file.\n";
unlink "$binDir/STOP";
print "All done.\n" . $stats->Show();


sub StopFile {
    open (my $oh, '>', "$binDir/STOP") || die "Could not open stop file: $!";
    print $oh "STOP\n";
    close $oh;
}

sub StartJob {
    my ($binDir, $dir, $extras) = @_;
    my $nohup = ($FIG_Config::win_mode ? "" : "nohup ");
    my $subDir = "$binDir/$dir";
    my $cmd = "bins_sample_pipeline $extras $dir $subDir >$subDir/run.log 2>$subDir/err.log";
    my $rc = system("$nohup$cmd &");
    my $time = scalar(localtime);
    # Check the start marker.
    my $start = 'Started';
    if (! -f "$subDir/contigs.fasta") {
        $start = 'Initiated';
    } elsif (-f "$subDir/START") {
        $start = 'Restarted';
    } else {
        # Create the start marker.
        open(my $oh, '>', "$subDir/START") || die "Could not create start marker for $subDir: $!";
        print $oh "\n";
    }
    print "$start job for $dir at $time.\n";
    $stats->Add("jobs$start" => 1);
}


