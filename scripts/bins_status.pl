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
use Bin::Package;
use FIG_Config;

=head1 Check Status of Bin Pipeline

    bins_status.pl [ options ] directory

Check the subdirectories of a specified directory to determine the progress of the binning runs. The status possibilities
are

=over 4

=item Downloaded

The reads have been downloaded. If the project type is C<NCBI>, this means C<site.tbl> exists.

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

Expectation processing completed-- C<expect.report.txt> exists or C<bins.rast.json> exists and there is no
abundance file.

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

=item terse

If specified, completed samples will not be shown, only samples in progress.

=item rerun

Used when rerunning assemblies.  Does not display assembled bins that are not running.

=item fix

If c<all>, all unprocessed samples will be removed. If C<empty>, all empty samples will be removed. The default is C<none>.

=item backout

Back out incomplete assemblies. This means removing the C<Assembly> directory so that the assembly can be restarted.

=item engine

Type of binning engine to use-- C<s> for the standard binner, C<2> for the alternate binner.

=item noIndex

If specified, the annotated genomes will not be indexed in PATRIC.

=back

=cut

# Files to keep during a reset.
use constant KEEPERS => { 'site.tbl' => 1, 'run.log' => 1, 'err.log' => 1, 'exclude.tbl' => 1, 'contigs.fasta' => 1,
        'output.contigs2reads.txt' => 1 };

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory',
                ['clean', 'clean up assembly information for complete samples'],
                ['resume', 'restart failed jobs'],
                ['terse', 'do not show completed bins'],
                ['rerun', 'do not show assembled bins that are not running'],
                ['project=s', 'project type for binning jobs', { default => 'MH' }],
                ['fix=s', 'remove unprocessed sample directories', { default => 'none' }],
                ['backout', 'back out incomplete assemblies'],
                ['maxResume=i', 'maximum number of running jobs for resume', { default => 20 }],
                ['engine=s', 'type of binning engine to use', { default => 's' }],
                ['noIndex', 'do not index bins in PATRIC'],
                ['reset', 'delete all binning results to force re-binning of all directories'],
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
my $proj = $opt->project;
my $fix = $opt->fix;
my $engine = $opt->engine;
my $noIndex = ($opt->noindex ? '--noIndex ' : '');
my $resetOpt = $opt->reset;
# Get a hash of the running subdirectories (Unix only).
my %running;
if (! $FIG_Config::win_mode) {
    my @jobs = `ps -AF`;
    for my $job (@jobs) {
        if ($job =~ /bins_sample_pipeline\s+(?:--\S+\s+)*(\w+)/) {
            $running{$1} = 1;
        }
    }
}
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
# Get the number of jobs we can resume. Each running job decrements this.
my $resumeLeft = $opt->maxresume - scalar (keys %running);
print "$resumeLeft job slots available.\n";
# Groups for printing.
my (@done, @downloaded, @other);
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    my $cleaned = (-d "$subDir/Assembly" ? "" : "  Assembly cleaned.");
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
            $stats->Add("site-$2" => 1);
        } else {
            $site = "Invalid";
        }
    } else {
        $site = "Error";
    }
    my $run = '';
    if ($running{$dir}) {
        $run = ', running';
    }
    my $label = "$subDir ($site$run)";
    # Are we resetting?
    if ($resetOpt && ! $run) {
        # Yes. Get the list of files and delete the binning stuff.
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
        $cleaned .= "  $count of $total files deleted by reset.";
    }
    # Check for the RAST completion file.
    my $rastFound = (-f "$subDir/bins.report.txt");
    # Check for the evaluation.
    my $evalDone = (-s "$subDir/Eval/index.tbl");
    # Determine the status.
    if (-s "$subDir/expect.report.txt") {
        $done = "Expectations Computed.";
    } elsif ($evalDone && ! -s "$subDir/$dir" . '_abundance_table.tsv') {
        $done = "Done (No Expectations).";
    } elsif ($evalDone) {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } else {
            push @other, "$label: Eval Complete.\n";
            $stats->Add(dirs7RastComplete => 1);
        }
    } elsif ($rastFound) {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } else {
            push @other, "$label: RAST Complete.\n";
            $stats->Add(dirs6RastComplete => 1);
        }
    } elsif (-s "$subDir/bin1.gto") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } else {
            my $i = 2;
            while (-s "$subDir/bin$i.gto") { $i++; }
            my $bins = $i - 1;
            push @other, "$label: RAST in Progress. $bins completed.\n";
            $stats->Add(binsAccumulating => $bins);
        }
        $stats->Add(dirs5RastPartial => 1);
    } elsif (-s "$subDir/bins.json") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } else {
            push @other, "$label: Bins Computed.\n";
        }
        $stats->Add(dirs4Binned => 1);
    } elsif (-f "$subDir/bins.report.txt") {
        $stats->Add(noBinsFound => 1);
        $done = "No bins found.";
    } elsif (-s "$subDir/sample.fasta") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } else {
            push @other, "$label: Binning in Progress.\n";
        }
        $stats->Add(dirs3Binning => 1);
    } elsif (-s "$subDir/contigs.fasta") {
        if (! $run && $opt->resume && $resumeLeft) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
            $resumeLeft--;
        } elsif (! $run && ! $opt->rerun) {
            push @other, "$label: Assembled.\n";
        }
        $stats->Add(dirs2Assembled => 1);
    } elsif (-d "$subDir/Assembly") {
        if (! $run && $opt->backout) {
            File::Copy::Recursive::pathrmdir("$subDir/Assembly");
            $stats->Add(assemblyBackout => 1);
            push @other, "$label: Downloaded.\n";
            $stats->Add(dirs0Downloaded => 1);
        } else {
            push @other, "$label: Assembling.\n";
            $stats->Add(dirs1Assembling => 1);
        }
    } elsif ($proj eq 'NCBI' && ! $sited) {
        push @other, "$label: Downloading.\n";
        $stats->Add('dirs.Downloading' => 1);
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
            StartJob($dir, $subDir, $gz, 'Started', $label, $proj);
            $runCount--;
        } elsif (! $run && $fix eq 'all') {
            # It's valid, and we are deleting it.
            File::Copy::Recursive::pathrmdir($subDir);
            $stats->Add(dirsDeleted => 1);
        } else {
            # It's valid, but we are leaving it alone.
            push @downloaded, "$label: Downloaded.\n";
            $stats->Add(dirs0Downloaded => 1);
        }
    }
    # If we are done, we process here and check for cleaning.
    if ($done) {
        my $show = ($opt->terse ? 0 : 1);
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
        # Check for evaluation results.
        if (open(my $ih, "$subDir/Eval/index.tbl")) {
            my $line = <$ih>;
            my ($good, $tot) = (0, 0);
            while (! eof $ih) {
                $line = <$ih>;
                if ($line =~ /\t1$/) {
                    $good++;
                    $stats->Add(goodBin => 1);
                }
                $tot++;
                $stats->Add(totBin => 1);
            }
            $cleaned .= "  $good of $tot bins.";
        }
        if ($show) {
            push @done, "$label: $done$cleaned\n";
        }
    }
}
print @done, @downloaded, @other;
print "\nAll done:\n" . $stats->Show();


sub StartJob {
    my ($dir, $subDir, $gz, $start, $label, $proj) = @_;
    my $realProj = ($proj eq 'NCBI' ? 'MH' : $proj);
    my $cmd = "bins_sample_pipeline --project=$realProj --engine=$engine $noIndex $gz $dir $subDir >$subDir/run.log 2>$subDir/err.log";
    my $rc = system("nohup $cmd &");
    push @other, "$label: $start $cmd.\n";
    $stats->Add("dirs0$start" => 1);
}