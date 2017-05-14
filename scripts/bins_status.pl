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

Project type for binning jobs. (C<MH>, C<AG>, C<HMP>, C<SYNTH>, or C<CT>)

=item run

Specifies that directories in a Downloaded status should be binned. The value should be a number. That number of
binning pipelines will be started; all subsequent directories will be left in a Downloaded state.

=item resume

Restart non-running jobs in progress.

=item terse

If specified, completed samples will not be shown, only samples in progress.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('directory',
                ['clean', 'clean up assembly information for complete samples'],
                ['resume', 'restart failed jobs'],
                ['terse', 'do not show completed bins'],
                ['project=s', 'project type for binning jobs', { required => 1 }],
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
my $runCount = $opt->run;
my $proj = $opt->project;
# Get a hash of the running subdirectories.
my %running;
my @jobs = `ps -AF`;
for my $job (@jobs) {
    if ($job =~ /bins_sample_pipeline\s+(?:--\S+\s+)*(\w+)/) {
        $running{$1} = 1;
    }
}
# Loop through the subdirectories.
opendir(my $ih, $directory) || die "Could not open directory $directory.";
my @dirs = sort grep { substr($_,0,1) ne '.' && -d "$directory/$_" } readdir $ih;
print scalar(@dirs) . " subdirectories found.\n";
# Groups for printing.
my (@done, @downloaded, @other);
for my $dir (@dirs) {
    $stats->Add(dirsTotal => 1);
    my $subDir = "$directory/$dir";
    my $rastFound = (-s "$subDir/bins.rast.json");
    my $cleaned = (-d "$subDir/Assembly" ? "" : "  Assembly cleaned.");
    my $done;
    # Determine the site.
    my $site;
    if (! -s "$subDir/site.tbl") {
        $site = "Unspecified";
    } elsif (open(my $sh, '<', "$subDir/site.tbl")) {
        my $line = <$sh>;
        if ($line =~ /(\S+)\t([^\t]+)\t(.+)/) {
            $site = "$1 $3";
            $stats->Add("site-$2" => 1);
        } else {
            $site = "Invalid";
        }
    } else {
        $site = "Error";
    }
    my $run = ($running{$dir} ? ', running' : '');
    my $label = "$subDir ($site$run)";
    # Determine the status.
    if (-s "$subDir/expect.report.txt") {
        $done = "Expectations Computed.";
    } elsif ($rastFound && ! -s "$subDir/$dir" . '_abundance_table.tsv') {
        $done = "Done (No Expectations).";
    } elsif ($rastFound) {
        if (! $run && $opt->resume) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
        } else {
            push @other, "$label: RAST Complete.\n";
            $stats->Add(dirs6RastComplete => 1);
        }
    } elsif (-s "$subDir/bin1.gto") {
        if (! $run && $opt->resume) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
        } else {
            my $i = 2;
            while (-s "$subDir/bin$i.gto") { $i++; }
            my $bins = $i - 1;
            push @other, "$label: RAST in Progress. $bins completed.\n";
        }
        $stats->Add(dirs5RastPartial => 1);
    } elsif (-s "$subDir/bins.json") {
        if (! $run && $opt->resume) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
        } else {
            push @other, "$label: Bins Computed.\n";
        }
        $stats->Add(dirs4Binned => 1);
    } elsif (-s "$subDir/bins.report.txt") {
        $stats->Add(noBinsFound => 1);
        $done = "No bins found.";
    } elsif (-s "$subDir/sample.fasta") {
        if (! $run && $opt->resume) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
        } else {
            push @other, "$label: Binning in Progress.\n";
        }
        $stats->Add(dirs3Binning => 1);
    } elsif (-s "$subDir/contigs.fasta") {
        if (! $run && $opt->resume) {
            StartJob($dir, $subDir, '', 'Restarted', $label, $proj);
        } else {
            push @other, "$label: Assembled.\n";
        }
        $stats->Add(dirs2Assembled => 1);
    } elsif (-d "$subDir/Assembly") {
        push @other, "$label: Assembling.\n";
        $stats->Add(dirs1Assembling => 1);
    } else {
        # Here the directory is downloaded. We may need to run the pipeline.
        my $startString = "";
        if ($runCount) {
            # Check for gz files.
            opendir(my $dh, $subDir) || die "Could not open directory $subDir: $!";
            my $found = grep { $_ =~ /\.fastq\.gz$/ } readdir $dh;
            closedir $dh;
            my $gz = ($found ? '--gz' : '');
            StartJob($dir, $subDir, $gz, 'Started', $label, $proj);
            $runCount--;
        } else {
            push @downloaded, "$label: Downloaded.\n";
            $stats->Add(dirs0Downloaded => 1);
        }
    }
    # If we are done, we process here and check for cleaning.
    if ($done) {
        my $show = ($opt->terse ? 0 : 1);
        $stats->Add(dirs7Done => 1);
        if ($clean && ! $cleaned) {
            print "Cleaning $subDir.\n";
            Bin::Package::CreateFromSample($subDir, $dir, $stats, 0, "$FIG_Config::data/GenomePackages");
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
            push @done, "$label: $done$cleaned\n";
        }
    }
}
print @done, @downloaded, @other;
print "\nAll done:\n" . $stats->Show();


sub StartJob {
    my ($dir, $subDir, $gz, $start, $label, $proj) = @_;
    my $cmd = "bins_sample_pipeline --project=$proj $gz $dir $subDir >$subDir/run.log 2>$subDir/err.log";
    my $rc = system("nohup $cmd &");
    push @other, "$label: $start $cmd.\n";
    $stats->Add("dirs0$start" => 1);
}