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
use Math::Round;
use SamplePipeline;

=head1 Assembly All Samples in a Binning Directory

    bins_assemble.pl [ options ] binDir

This is a long-running job that assembles jobs for binning by L<bins_manage.pl>.  It will look for one assembly at
a time.  The assembly cannot have been started, which means there cannot be an C<Assembly> directory.  If all
samples are assembled it will stop.  It can also be stopped using the C<--stopAsm> option on L<bins_status.pl>.

=head2 Parameters

The positional parameter is the name of the binning directory.  All the samples must be in subdirectories.

The following command-line options are supported.

=over 4

=item small

Only assemble small samples.

=back

=cut


$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
        ['small', 'only assemble small samples']
        );
# Get the input directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No binning directory specified.";
} elsif (! -d $binDir) {
    die "Binning directory $binDir missing or invalid.";
}
# Here we want to go into the loop.
# Loop until we find the stop file.
while (! -f "$binDir/STOPASM") {
    # This will count the incomplete directories.
    my $incomplete = 0;
    # Get the subdirectories.  We filter out directories with contig files, but we don't filter on assemblies in
    # progress, because we want to do that as close as possible to the point where the assembly is started.
    opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
    my @samples = sort grep { substr($_,0,1) ne '.' && -d "$binDir/$_" && ! -s "$binDir/$_/contigs.fasta" } readdir $dh;
    closedir $dh;
    # This will be set to TRUE if we found a directory.
    my $found;
    # Loop through the directories.
    while (! $found && scalar @samples) {
        my $sample = pop @samples;
        # Can we assemble this directory?
        my $subDir = "$binDir/$sample";
        if (! -d "$subDir/Assembly") {
            File::Copy::Recursive::pathmk("$subDir/Assembly") || die "Could not claim $subDir: $!";
            # We have claimed this directory.  Mark it as having an assembly in process.
            open(my $oh, '>', "$subDir/ASSEMBLE") || die "Could not create assembly marker for $sample: $!";
            print $oh "\n";
            close $oh;
            # Create the options hash for the assembly pipeline call.
            my %options = (noBin => 1);
            SamplePipeline::PrepareAssembly($subDir, \%options);
            if ($opt->small && $options{large}) {
                print "Skipping $sample:  too big for this machine.\n";
                # Release the directory.
                File::Copy::Recursive::pathempty("$subDir/Assembly") || die "Could not clean up $subDir: $!";
                rmdir "$subDir/Assembly";
            } else {
                # Perform the assembly.
                SamplePipeline::Process($subDir, %options);
                # Erase the FASTQ files if we were successful.
                if (-s "$subDir/contigs.fasta") {
                    SamplePipeline::ClearAssembly($subDir, $sample);
                }
            }
            # Remove the marker.
            unlink "$subDir/ASSEMBLE";
            # Denote we found something to assemble.
            $found = 1;
        }
    }
    # Did we assemble something?
    if (! $found) {
        # No.  Stop the loop.
        open(my $oh, '>', "$binDir/STOPASM") || die "Could not create stop file: $!";
        print $oh "\n";
    }
}