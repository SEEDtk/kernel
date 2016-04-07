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
use Bin::Compute;
use Bin::Score;
use Bin::Analyze;
use Shrub;
use Stats;
use Time::HiRes qw(time);

=head1 Build Bins From Community Contigs Using Multiple Scoring Schemes (Algorithm 2)

    bins_decide [ options ] workDirectory <contigFile

This script reads contig information from the standard input and partitions them into bins based
on criteria relating to closest reference genomes, universal roles, coverage, and tetranucleotides.
The contigs are represented by L<Bin> objects.

Unlike L<bins_compute.pl>, which uses a single weighted scoring scheme, this method takes a list of
scoring schemes as input and clusters contigs based on what portion of the schemes predict they should
go together. It therefore runs the standard binning algorithm many times and tallies how often two
contigs are put together.

=head2 Parameters

The positional parameters are the name of a working directory to contain temporary and output files and the
name of a scheme file (see below).

The command-line options are those found in L<ScriptUtils/ih_options> and L<Shrub/script_options> plus the following.

=over 4

=item cutoff

Each scheme produces a YES vote or a NO vote for putting two contigs into the same bin. This parameter is the percentage of
YES votes that are required to put two contigs in the same bin during the final assembly.

=item scoreFile

If specified, is presumed to be a file containing connection counts, tab-delimited, with each record consisting of two
contig IDs followed by a count. In the presence of this file, the initial binning process can be skipped, and the
output produced simply by comparing the counts to the cutoff.

=back

=head2 Input Files

The main input file contains one or more L<Bin> objects sequentially in the format described by L<Bin/Bin Exchange Format>.
Each represents a single contig from the input community.

The scheme file contains one or more records representing schemes for scoring pairs of contigs. Each record is tab-delimited,
and consists of (0) a rank, (1) a label, (2) the coverage weight, (3) the tetranucleotide weight, (4) the closest-reference-genome
weight, (5) the penalty for duplicate universal roles, (6) the universal role weight, and (7) the unscaled minimum score.
The scheme file is designed to be taken from the output of L<ga.pl>.

=head3 Output Files

The output files are as follows, all in the working directory.

=over 4

=item bins.json

A file of L<Bin> objects in JSON format, suitable for reading by L<Bin/ReadBins>. These represent the computed bins.

=item connections.tbl

A file of connection counts for each pair of contigs. The file is tab-delimited, each line consisting of two contig IDs followed
by the count.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDirectory schemeFile', ScriptUtils::ih_options(), Shrub::script_options(),
        ['cutoff|c=f', 'percent of schemes that must approve connecting two contigs', { default => 70 }],
        ['scoreFile|s=s', 'connection score inptu file']
        );
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the list of contigs. These are read in as bins.
print "Reading contigs from input.\n";
my $binList = Bin::ReadContigs($ih);
# Verify the working directory and scheme file.
my ($workDir, $schemeFile) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.";
}
if (! $schemeFile) {
    die "No scheme file specified.";
} elsif (! -f $schemeFile) {
    die "Scheme file $schemeFile not found.";
}
# Create an empty scoring object.
my $score = Bin::Score->new(0, 0, 0, 0, 0, 0);
# Create the computation object.
my $computer = Bin::Compute->new($score);
# This will count the number of schemes.
my $schemes = 0;
# This will contain the final scoring vector.
my @scores;
# Open the scheme file.
open(my $sh, '<', $schemeFile) || die "Could not open $schemeFile: $!";
# Do we have a score file?
if ($opt->scorefile) {
    open(my $kh, '<', $opt->scorefile) || die "Could not open scoreFile: $!";
    # Count the schemes.
    while (! eof $sh) {
        my $line = <$sh>;
        $schemes++;
    }
    # Compute the actual numerical cutoff. Any score greater than or equal to this is kept.
    my $thresh = int(($opt->cutoff * $schemes + 99) / 100);
    # Read the scores into the score vector.
    while (! eof $kh) {
        my $line = <$kh>;
        chomp $line;
        my ($contig1, $contig2, $count) = split /\t/, $line;
        if ($count >= $thresh) {
            push @scores, [$contig1, $contig2, $count];
            $stats->Add(scoreKept => 1);
        } else {
            $stats->Add(scoreDiscarded => 1);
        }
    }
} else {
    # This will contain the statistics generated by the subprograms.
    my $subStats;
    # For each contig, this hash will map the bin ID to a hash counting the number of times each other contig was predicted
    # to belong to it. The first contig will always be lexically earlier. Thus there would be a value for $binTogether{ABC}{BBC}
    # but not $binTogether{BBC}{ABC}.
    my %binTogether;
    # Loop through the schemes.
    while (! eof $sh) {
        # Read the scheme.
        my $line = <$sh>;
        chomp $line;
        my (undef, $label, $covgweight, $tetraweight, $refweight, $unipenalty, $uniweight, $inMinscore) = split /\t/, $line;
        my $minscore = Bin::Score::scale_min($covgweight, $tetraweight, $refweight, $uniweight, $inMinscore);
        $schemes++;
        print "Processing scheme #$schemes: $label (covg = $covgweight, tetra = $tetraweight, ref = $refweight, uni = $uniweight penalty $unipenalty, min = $minscore).\n";
        # Store it in the scoring object.
        $score->Reset($covgweight, $tetraweight, $refweight, $unipenalty, $uniweight, $minscore);
        # Compute the bins.
        my $start = time();
        my @startBins = map { Bin->new_copy($_) } @$binList;
        my $bins = $computer->ProcessScores(\@startBins);
        my $duration = time() - $start;
        $stats->Accumulate($subStats);
        $stats->Add(duration => int($duration + 0.5));
        # Count the connections.
        print "Counting connections.\n";
        for my $bin (@$bins) {
            # Get the sorted list of contigs. Each contig has a connection to all the contigs after it (because of the sorting).
            my @contigs = sort $bin->contigs;
            # We shift off the first contig, then connect it to all the remaining ones. This process continues until all the contigs
            # are processed. Note if there is only one contig left, there is nothing to compare it to.
            while (scalar(@contigs) > 1) {
                my $contig1 = shift @contigs;
                for my $contig2 (@contigs) {
                    $binTogether{$contig1}{$contig2}++;
                    $stats->Add(connectionsCounted => 1);
                }
            }
        }
    }
    # Compute the actual numerical cutoff. Any count greater than or equal to this is kept.
    my $thresh = int(($opt->cutoff * $schemes + 99) / 100);
    # Convert the contig hash into a score vector and checkpoint it to the file.
    my $scoreFile = "$workDir/connections.tbl";
    open(my $sch, '>', $scoreFile) || die "Could not open $scoreFile: $!";
    print "Creating score vector.\n";
    for my $contig1 (keys %binTogether) {
        my $contigH = $binTogether{$contig1};
        for my $contig2 (keys %$contigH) {
            my $count = $contigH->{$contig2};
            print $sch "$contig1\t$contig2\t$count\n";
            if ($count >= $thresh) {
                $stats->Add(connectionKept => 1);
                push @scores, [$contig1, $contig2, $count];
            } else {
                $stats->Add(connectionDiscarded => 1);
            }
        }
    }
    close $sch;

}
# Cluster the contigs based on the scores.
my $bins = $computer->ClusterContigs($binList, \@scores, closed => 1);
# Write the resulting bins.
print "Writing bins.\n";
open(my $oh, ">", "$workDir/bins.json") || die "Could not open bins output file: $!";
for my $bin (@$bins) {
    $bin->Write($oh);
}
close $oh;
# Analyze the bins.
my $quality = Bin::Analyze::Quality($shrub, $bins);
print "Quality score = $quality.\n";
my $report = Bin::Analyze::Report($shrub, $bins);
print "Quality Report\n" . $report->Show() . "\n";
# Output the statistics.
print "All done.\n" . $stats->Show();
