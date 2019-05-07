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

This script will look at an evaluated binning directory and try to find the representatives of the good genomes found.
A report on the genomes found will be placed in the C<reps.tbl> file in the binning directory.

=head2 Parameters

The positional parameters are the name of the binning directory and the name of the representative-genome directory
(described in L<RepGenomeDb>).

Additional command-line options are the following.

=over 4

=item recursive

If specified, the specified input directory is presumed to contain multiple subdirectories that are binning directories.

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
    ['recursive|R', 'process binning subdirectories']);
# Create the statistics object.
my $stats = Stats->new();
# This will be the option hash for creating the GEOs.
my %options = (stats => $stats, detail => 0, logH => \*STDOUT);
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
print "Loading representative genome data from $repDir.\n";
my $repDB = RepGenomeDb->new_from_dir($repDir, unconnected => 1, verbose => 1);
# Now we need to find the binning data.
my @binDirs = ($binDir);
if ($opt->recursive) {
    # Here we are reading subdirectories.
    print "Scanning $binDir for completed binning directories.\n";
    opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
    @binDirs = map { "$binDir/$_" } grep { substr($_, 0, 1) ne '.' && -s "$binDir/$_/Eval/index.tbl" } readdir $dh;
    closedir $dh;
    print scalar(@binDirs) . " completed binning directories found.\n";
}
# Get the cutoff score for this RepDB.
my $minScore = $repDB->score();
# Loop through the binning directories.
for my $dir (@binDirs) {
    # Get the GEOs for this bin.
    my $geos = GeoGroup->new(\%options, $dir);
    $stats->Add(dirsIn => 1);
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
                    print "$genome $name is an outlier.\n";
                    $stats->Add(outlier => 1);
                    $class = 'outlier';
                } else {
                    $stats->Add(represented => 1);
                    $class = '';
                }
            }
            P3Utils::print_cols([$genome, $name, $repID, $rep_name, $score, $class], oh => $oh);
        }
    }
}
print "All done.\n" . $stats->Show();

