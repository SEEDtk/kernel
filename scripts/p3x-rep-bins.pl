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

=head1 Check for Unrepresented Genomes in Bins

    p3x-rep-bins.pl [options] binDir repDir

This script will look at evaluated binning directories and try to find the representatives of the good genomes found.
A report on the genomes found will be placed in the C<reps.tbl> file in the binning directory and in the standard
output.

=head2 Parameters

The positional parameters are the name of the master binning directory (containing the actual binning directory as
subdirectories) and the name of the representative-genome directory (described in L<RepGenomeDb>).

Additional command-line options are the following.

=over 4

=item missing

If specified, directories with an existing C<reps.tbl> file will not be reprocessed.

=item verbose

Display status messages on STDERR.

=back

=cut

use strict;
use GEO;
use P3Utils;
use RepGenomeDb;
use GeoGroup;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('binDir repDir',
    ['missing', 'only process new directories'],
    ['verbose|debug|v', 'display progress messages on STDERR']);
# Create the statistics object.
my $stats = Stats->new();
# Get the options.
my $debug = $opt->verbose;
my $missing = $opt->missing;
# This will be the option hash for creating the GEOs.
my %options = (stats => $stats, detail => 0, logH => \*STDERR);
# Get the positional parameters.
my ($binDir, $repDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "Input directory $binDir invalid or missing.";
} elsif (! $repDir) {
    die "No rep-genomes directory specified.";
} elsif (! -d $repDir) {
    die "Rep-genomes directory $repDir invalid or missing.";
} elsif (! -s "$repDir/6.1.1.20.fasta") {
    die "$repDir does not appear to be a representative-genomes directory.";
}
# Create the rep-genomes object.
print STDERR "Loading representative genome data from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir, unconnected => 1, verbose => 1);
# Now we need to find the binning data.
# Here we are reading subdirectories.
print "Scanning $binDir for completed binning directories.\n";
opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
my @binDirs = grep { substr($_, 0, 1) ne '.' && -s "$binDir/$_/Eval/index.tbl" } readdir $dh;
closedir $dh;
print STDERR scalar(@binDirs) . " completed binning directories found.\n" if $debug;
# Get the cutoff score for this RepDB.
my $minScore = $repDB->score();
print STDERR "Representation score is $minScore.\n" if $debug;
# Start the main output file.
P3Utils::print_cols(['directory', 'genome', 'name', 'rep', 'rep_name', 'score', 'class']);
# Loop through the binning directories.
for my $sample (@binDirs) {
    $stats->Add(dirsIn => 1);
    my $dir = "$binDir/$sample";
    if ($missing && -s "$dir/reps.tbl") {
        # Here we do not need to re-process the directory, but we need to include its data in the main
        # report and the stats.
        print STDERR "Reading report for $sample.\n";
        open(my $ih, '<', "$dir/reps.tbl") || die "Could not open reps.tbl for $dir: $!";
        # Skip the header.
        my $line = <$ih>;
        # Process the data lines.
        while (! eof $ih) {
            my $line = <$ih>;
            print "$sample\t$line";
            if ($line =~ /extreme$/) {
                $stats->Add(extremeOutlier => 1);
            } elsif ($line =~ /outlier$/) {
                $stats->Add(outlier => 1);
            } else {
                $stats->Add(represented => 1);
            }
        }
    } else {
        # Get the GEOs for this bin.
        my $geos = GeoGroup->new(\%options, $dir);
        # Create the output file.
        open(my $oh, '>', "$dir/reps.tbl") || die "Could not open output for $binDir: $!";
        P3Utils::print_cols(['genome', 'name', 'rep', 'rep_name', 'score', 'class'], oh => $oh);
        # Loop through the GEOs.
        my $geoList = $geos->geoList;
        for my $geo (@$geoList) {
            # Only process good genomes.
            if (! $geo->is_good) {
                $stats->Add(badGenome => 1);
            } else {
                $stats->Add(goodGenome => 1);
                # These will be the output values.
                my ($genome, $name, $repID, $rep_name, $score, $class) = ($geo->id, $geo->name, '', '', 0, 'extreme');
                ($repID, $score) = $repDB->find_rep($geo->seed);
                if (! $repID) {
                    print "$genome $name is an extreme outlier.\n";
                    $stats->Add(extremeOutlier => 1);
                } else {
                    $rep_name = $repDB->rep_object($repID)->name;
                    if ($score < $minScore) {
                        print "$genome $name is an outlier with score $score.\n";
                        $stats->Add(outlier => 1);
                        $class = 'outlier';
                    } else {
                        $stats->Add(represented => 1);
                        $class = '';
                    }
                }
                P3Utils::print_cols([$genome, $name, $repID, $rep_name, $score, $class], oh => $oh);
                P3Utils::print_cols([$sample, $genome, $name, $repID, $rep_name, $score, $class]);
            }
        }
    }
}
print STDERR "All done.\n" . $stats->Show();

