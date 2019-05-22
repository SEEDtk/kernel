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
use FIG_Config;
use Math::Round;

=head1 Check Status of Bin Pipeline

    bins_status.pl [ options ] directory

Check the subdirectories of a specified directory to determine the progress of the binning runs. The status possibilities
are

=over 4

=item Downloaded

The reads have been downloaded. This means C<site.tbl> exists.

=item Assembling

The reads are being assembled into contigs-- C<Assembly> exists.

=item Assembled

The reads have been assembled into contigs-- C<Assembly/contigs.fasta> exists.

=item Bins Computed

The contigs have been sorted into bins-- C<bins.json> exists.

=item RAST in Progress

One or more bins have been processed by RAST-- C<bin>I<X>C<.gto> exists for one or more bin numbers I<X>.

=item RAST Complete

RAST processing completed-- C<bins.rast.json> exists.

=item Done

Evaluation processing completed-- C<Eval/index.tbl> exists.

=back

=head2 Parameters

The single positional parameter is the name of the directory containing the sample sub-directories.

The following command-line options are supported.

=over 4

=item clean

Remove assembly information for completed jobs.

=item project

Project type for binning jobs. (C<MH>, C<AG>, C<HMP>, C<SYNTH>, C<NCBI>, or C<CT>)

=item run

Specifies that directories in a Downloaded status should be binned. The value should be a number. That number of
binning pipelines will be started; all subsequent directories will be left in a Downloaded state.

=item resume

Restart non-running jobs in progress.

=item maxResume

Maximum number of running jobs to allow when resuming. The default is C<20>. This is always approximate: it is just a
stopgap to avoid overwhelming the machine.

=item all

If specified, completed samples, and assembled or downloaded bins that are not running will be shown.

=item fix

If c<all>, all unprocessed samples will be removed. If C<empty>, all empty samples will be removed. The default is C<none>.

=item backout

Back out incomplete assemblies. This means removing the C<Assembly> directory so that the assembly can be restarted.

=item engine

Type of binning engine to use-- C<s> for the standard binner, C<2> for the alternate binner.

=item noIndex

If specified, the annotated genomes will not be indexed in PATRIC.

=item rebin

Remove binning results from samples stopped during evaluation.  If the value is C<all>, also remove binning results from
samples stopped during binning.

=item stopFile

The name of a file to contain error information about stopped jobs.  The default is C<stoppedJobs.log> in the
current directory.

=item target

The target number of assemblies to use for predicting the next command.  The default is C<4>.

=item stop

If specified, a file will be generated to stop the L<bins_manage.pl> process.

=back

=cut

# Files to keep during a reset.
use constant KEEPERS => { 'site.tbl' => 1, 'run.log' => 1, 'err.log' => 1, 'exclude.tbl' => 1, 'contigs.fasta' => 1,
        'output.contigs2reads.txt' => 1 };

# Save the parameters.
my @saved = grep { ! ($_ =~ /^--(?:clean|resume|fix|backout|reset|run)/) } @ARGV;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory',
                ['clean', 'clean up assembly information for complete samples'],
                ['resume', 'restart failed jobs'],
                ['all', 'show all directories'],
                ['fix=s', 'remove unprocessed sample directories', { default => 'none' }],
                ['backout', 'back out incomplete assemblies'],
                ['maxResume=i', 'maximum number of running jobs for resume', { default => 20 }],
                ['engine=s', 'type of binning engine to use', { default => 's' }],
                ['noIndex', 'do not index bins in PATRIC'],
                ['reset', 'delete all binning results to force re-binning of all directories'],
                ['rebin:s', 'reset samples that are stopped during evaluation and/or binning'],
                ['stopFile=s', 'file to contain stopped-job error data', { default => 'stoppedJobs.log' }],
                ['stop', 'stop the binning manager loop'],
                ['stopAsm', 'stop the binning assembler loop'],
                ['target=i', 'maximum number of assembly jobs', { default => 4 }],
                ['run=i', 'run binning pipeline on new directories', { default => 0 }]);
my $stats = Stats->new();
# Get the main directory name.
my ($directory) = @ARGV;
if (! $directory) {
    die "No directory specified.";
} elsif (! -d $directory) {
    die "Invalid directory name $directory.";
}
# Save the options.
my $clean = $opt->clean;
my $runCount = $opt->run // 0;
my $fix = $opt->fix;
my $engine = $opt->engine;
my $noIndex = ($opt->noindex ? '--noIndex ' : '');
my $resetOpt = $opt->reset;
# Make sure the rebin value is comparable and easy to discern.
my $rebin = $opt->rebin;
if (! defined $rebin) {
    $rebin = '';
} elsif (! $rebin) {
    $rebin = 'some';
}
# Process the manager-stop options first.
if ($opt->stop) {
    open (my $oh, '>', "$directory/STOP") || die "Could not open stop file: $!";
    print $oh "STOP\n";
    close $oh;
    print "STOP file created.\n";
}
if ($opt->stopasm) {
    open (my $oh, '>', "$directory/STOPASM") || die "Could not open stop-assembly file: $!";
    print $oh "STOP\n";
    close $oh;
    print "STOPASM file created.\n";
}
# Get an output handle for the stopped-jobs file.
open(my $sh, '>', $opt->stopfile) || die "Could not open stopped-jobs file: $!";
# Get a hash of the running subdirectories (Unix only).
my %running;
if (! $FIG_Config::win_mode) {
    my @jobs = `ps -Af`;
    for my $job (@jobs) {
        if ($job =~ /bins_sample_pipeline\s+(?:--\S+\s+)*(\w+)/) {
            $running{$1} = 1;
        }
    }
}
# These will be used to determine production ratio and assembly total.
my ($totGood, $totDone, $asming, $stopped) = (0, 0, 0, 0);
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
# Get the number of jobs we can resume. Each running job decrements this.
my $resumeLeft = $opt->maxresume - scalar (keys %running);
print "$resumeLeft job slots available.\n";
# Subtract the number of jobs to start from the number of jobs left.
$resumeLeft -= $runCount;
# Groups for printing.
my (@done, @downloaded, @other);
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    my $cleaned = (-d "$subDir/Assembly" ? "" : "  Assembly cleaned.");
    my $binStatus = '';
    my $done;
    # Determine the site.
    my ($site, $sited);
    if (! -s "$subDir/site.tbl") {
        $site = "Unspecified";
    } elsif (open(my $sh, '<', "$subDir/site.tbl")) {
        my $line = <$sh>;
        if ($line =~ /(\S+)\t([^\t]+)\t(.+)/) {
            $site = "$1 $3";
            $sited = 1;
        } else {
            $site = "Invalid";
            $stats->Add(siteInvalid => 1);
        }
    } else {
        $site = "Error";
        $stats->Add(siteError => 1);
    }
    my $run = '';
    if ($running{$dir}) {
        $run = ', running';
    }
    my $label = "$subDir ($site$run)";
    # Are we resetting?
    if ($resetOpt && ! $run) {
        # Yes. Delete the binning stuff.
        ResetSample($dir, $subDir);
    }
    # Check for the RAST completion file.
    my $rastFound = (-f "$subDir/bins.report.txt");
    # Check for the evaluation.
    my $evalDone = (-s "$subDir/Eval/index.tbl");
    # Determine the status.
    if ($evalDone && ! -s "$subDir/bins.found.tbl") {
        $done = "Done (No bins found).";
        $stats->Add(noBinsFound => 1);
    } elsif ($evalDone && ! -s "$subDir/ref.genomes.scores.tbl") {
        $done = "Done (No binnable seeds).";
        $stats->Add(noBinsFound => 1);
    } elsif ($evalDone) {
        $done = "Done.";
        # Check for evaluation results.
        if (open(my $ih, "$subDir/Eval/index.tbl")) {
            my $line = <$ih>;
            my ($good, $tot) = (0, 0);
            while (! eof $ih) {
                $line = <$ih>;
                if ($line =~ /\t1$/) {
                    $good++;
                    $totGood++;
                    $stats->Add(binsGood => 1);
                } else {
                    $stats->Add(binsBad => 1);
                }
                $tot++;
                $stats->Add(binsTotal => 1);
            }
            $totDone++;
            if ($tot) {
                $binStatus = "  $good of $tot bins.";
            }
        }
    } elsif ($rastFound) {
        if (! $run && $rebin) {
            ResetSample($dir, $subDir);
        }
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label);
            $resumeLeft--;
        } else {
            push @other, "$label: RAST Complete.\n";
            $stats->Add(dirs7RastComplete => 1);
            if (! $run) {
                RecordStopped($sh, $label, 'Evaluation', $subDir, \$stopped);
            }
        }
    } elsif (-s "$subDir/bin1.gto") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label);
            $resumeLeft--;
        } else {
            my $i = 2;
            while (-s "$subDir/bin$i.gto") { $i++; }
            my $bins = $i - 1;
            push @other, "$label: RAST in Progress. $bins completed.\n";
            $stats->Add(binsAccumulating => $bins);
            RecordStopped($sh, $label, 'RAST Processing', $subDir, \$stopped) if (! $run);
        }
        $stats->Add(dirs6RastPartial => 1);
    } elsif (-s "$subDir/bins.json") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label);
            $resumeLeft--;
        } else {
            push @other, "$label: Bins Computed.\n";
            RecordStopped($sh, $label, 'RAST Processing', $subDir, \$stopped) if (! $run);
        }
        $stats->Add(dirs5Binned => 1);
    } elsif (-f "$subDir/bins.report.txt") {
        $stats->Add(noBinsFound => 1);
        $done = "No bins found.";
    } elsif (-s "$subDir/sample.fasta") {
        if (! $run && $rebin eq 'all') {
            # Stopped and we are re-binning these. Delete the binning stuff.
            ResetSample($dir, $subDir);
        }
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label);
            $resumeLeft--;
        } else {
            RecordStopped($sh, $label, 'Binning', $subDir, \$stopped) if (! $run);
            push @other, "$label: Binning in Progress.\n";
        }
        $stats->Add(dirs4Binning => 1);
    } elsif (-s "$subDir/contigs.fasta") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label);
            $resumeLeft--;
        } elsif (! $run && $opt->all) {
            push @other, "$label: Assembled.\n";
        }
        $stats->Add(dirs3Assembled => 1);
    } elsif (-d "$subDir/Assembly") {
        if (! $run && $opt->backout) {
            File::Copy::Recursive::pathrmdir("$subDir/Assembly");
            $stats->Add(assemblyBackout => 1);
            push @other, "$label: Downloaded.\n";
            $stats->Add(dirs1Downloaded => 1);
        } else {
            # Here we are assembling.  Get the time in progress if it is running or is on an assembly machine.
            my $duration = '';
            if ($run || -f "$subDir/ASSEMBLE") {
                $duration = -M "$subDir/Assembly/params.txt";
                if ($duration) {
                    $duration = "  " . Math::Round::nearest(0.1, 24 * $duration) . " hours.";
                } else {
                    $duration = '';
                }
                $asming++ if $run;
            }
            push @other, "$label: Assembling.$duration\n";
            $stats->Add(dirs2Assembling => 1);
        }
    } elsif (! -s "$subDir/site.tbl") {
        # A download is in progress here.  Compute the progress.
        opendir(my $dh, $subDir) || die "Could not open directory $subDir: $!";
        my $count = 0;
        map { $count += -s "$subDir/$_" } grep { $_ =~ /\.(?:fastq|fq)/ } readdir $dh;
        if ($count > 1000000) {
            $count = int($count/1000000) . "M";
        }
        push @other, "$label: Downloading. $count so far.\n";
        $stats->Add(dirs0Downloading => 1);
    } else {
        # Here the directory is downloaded. We may need to fix it or run the pipeline.
        opendir(my $dh, $subDir) || die "Could not open directory $subDir: $!";
        my @files = grep { $_ =~ /\.(?:fastq|fq)/ } readdir $dh;
        closedir $dh;
        my $found = scalar @files;
        if (! $found) {
            my $status = "$label: Empty.";
            $stats->Add(dirs0Empty => 1);
            if ($fix ne 'none') {
                File::Copy::Recursive::pathrmdir($subDir);
                $status .= "  Deleted.\n";
                $stats->Add(dirsDeleted => 1);
            } else {
                $status .= "\n";
            }
            push @other, $status;
        } elsif ($runCount > 0) {
            # It's valid, and we want to run it. Check for gz files.
            $found = grep { $_ =~ /q\.gz$/ } @files;
            my $gz = ($found ? '--gz' : '');
            StartJob($dir, $subDir, $gz, 'Started', $label);
            $runCount--;
        } elsif (! $run && $fix eq 'all') {
            # It's valid, and we are deleting it.
            File::Copy::Recursive::pathrmdir($subDir);
            $stats->Add(dirsDeleted => 1);
        } else {
            # It's valid, but we are leaving it alone.
            $stats->Add(dirs1Downloaded => 1);
            if ($opt->all) {
                push @downloaded, "$label: Downloaded.\n";
            }
        }
    }
    # If we are done, we process here and check for cleaning.
    if ($done) {
        my $show = $opt->all || 0;
        $stats->Add(dirs8Done => 1);
        if ($clean && ! $cleaned) {
            print "Cleaning $subDir.\n";
            $cleaned = "  Cleaning Assembly.";
            File::Copy::Recursive::pathempty("$subDir/Assembly");
            rmdir "$subDir/Assembly";
            opendir(my $dh, $subDir) || die "Could not open $subDir: $!";
            my @fastqs = grep { $_ =~ /\.(?:fastq|fq)$/ } readdir $dh;
            for my $fastq (@fastqs) {
                unlink "$subDir/$fastq";
            }
            $stats->Add(dirsCleaned => 1);
            $show = 1;
        }
        if ($show) {
            push @done, "$label: $done$cleaned$binStatus\n";
        }
    }
}
if ($totDone) {
    print "Production ratio is " . Math::Round::nearest(0.01, $totGood / $totDone) . ".\n";
}
if ($asming) {
    print "$asming assemblies in progress.\n";
} else {
    print "No assemblies in progress.\n";
}
if ($stopped) {
    print "$stopped jobs are stopped.\n";
}
print @done, @downloaded, @other;
print "\nAll done:\n" . $stats->Show();


sub StartJob {
    my ($dir, $subDir, $gz, $start, $label) = @_;
    my $cmd = "bins_sample_pipeline --engine=$engine $noIndex $gz $dir $subDir >$subDir/run.log 2>$subDir/err.log";
    my $rc = system("nohup $cmd &");
    push @other, "$label: $start $cmd.\n";
    $stats->Add("dirsX$start" => 1);
}

sub ResetSample {
    my ($dir, $subDir) = @_;
    my ($count, $total) = (0, 0);
    opendir(my $dh, $subDir) || die "Could not open work directory: $!";
    my @files = grep { -f "$subDir/$_" } readdir $dh;
    for my $file (@files) {
        my $fullName = "$subDir/$file";
        $total++;
        unless ($fullName =~ /_abundance_table.tsv$/ || $fullName =~ /\.fastq$/ || $fullName =~ /\.fq/ ||
                KEEPERS->{$file}) {
            unlink $fullName;
            $count++;
        }
    }
    if (-d "$subDir/Eval") {
        File::Copy::Recursive::pathrmdir("$subDir/Eval") || die "Could not remove Eval from $subDir: $!";
    }
    print "$count of $total files deleted by reset of $dir.\n";
}

sub RecordStopped {
    my ($sh, $label, $status, $subDir, $pStopped) = @_;
    print $sh "*** $label stopped during $status.\n";
    my $ih;
    my $size = -s "$subDir/err.log";
    if (! $size) {
        print $sh "  * Error log file is missing or empty.\n";
    } elsif (! open($ih, '<', "$subDir/err.log")) {
        print $sh "  * Could not open error log file: $!\n";
    } else {
        my $pos = $size - 500;
        if ($pos > 0) {
            seek $ih, $pos, 0;
            # Throw out the first line fragment.
            my $line = <$ih>;
        }
        # Echo the rest.
        while (! eof $ih) {
            my $line = <$ih>;
            print $sh "    $line";
        }
        print $sh "\n";
    }
    $$pStopped++;
}