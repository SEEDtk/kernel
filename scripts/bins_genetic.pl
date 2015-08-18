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
use Bin;
use Bin::Score;
use Bin::Compute;
use Bin::Analyze;
use AI::Genetic;
use Time::HiRes qw(time);
use IO::File;

=head1 Search for Good Binning Parameters

    bins_genetic.pl [ options ] workDirectory

This program runs the community binning process multiple times to attempt to find the best parameter values,
using the L<AI::Genetic> genetic algorithm.

=head2 Parameters

The single positional parameter is the name of a working directory to contain temporary and output files.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item generations

The number of generations for which to run the search.

=item population

The population size to use for the search.

=back

=head2 Input File

The input file contains one or more L<Bin> objects sequentially in the format described by L<Bin/Bin Exchange Format>.
Each represents a single contig from the input community.

=head3 Output Files

The output files are as follows, all in the working directory.

=over 4

=item genes.log

A text file containing the input scores and the analysis statitics of each run.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDirectory', ScriptUtils::ih_options(),
        ['generations', 'number of generations to run', { default => 50 }],
        ['population', 'starting population', { default => 100 }],
        );
# Create the genetic search object.
my $ga = new AI::Genetic(
        -fitness => \&ComputeBins,
        -type    => 'rangevector',
        -population => $opt->population,
        -crossover  => 0.9,
        -mutation   => 0.01,
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the list of contigs. These are read in as bins.
print "Reading contigs from input.\n";
my $binList = Bin::ReadContigs($ih);
print scalar(@$binList) . " contigs found.\n";
# Verify the working directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.\n";
}
# Open the log file. We take special pains so that it is unbuffered.
my $oh = IO::File->new(">$workDir/genes.log") || die "Could not open output log.";
$oh->autoflush(1);
# Initialize the scoring parameter ranges.
$ga->init([[0, 10], [1, 10], [0, 10], [0, 10], [0, 10], [5, 7]]);
# Search for the best scoring parameters.
for (my $g = 1; $g <= $opt->generations; $g++) {
    print "*** Generation $g.\n";
    my $start = time;
    # Process this generation.
    $ga->evolve('randomUniform', 1);
    my ($s) = $ga->getFittest(1);
    print "Best score is " . $s->score . ". Parms = [" . join(", ", @{$s->genes}) . "]. " . int(time - $start + 0.5) . " seconds/generation.\n";
}

# This is the fitness function. It takes the scoring parameters as input and produces a quality score.
sub ComputeBins {
    my ($scoreParms) = @_;
    my $start = time;
    # Scale the scores.
    my ($covgweight, $tetraweight, $refweight, $uniweight, $unipenalty, $minscore) = @$scoreParms;
    my $denom = ($covgweight + $tetraweight + $refweight + $uniweight);
    ($covgweight, $tetraweight, $refweight, $uniweight) = map { $_ / $denom } ($covgweight, $tetraweight, $refweight, $uniweight);
    $minscore /= 10;
    $unipenalty /= 5;
    my @realScores = ($covgweight, $tetraweight, $refweight, $unipenalty, $uniweight, $minscore);
    # Create the scoring object.
    my $score = Bin::Score->new(@realScores);
    # Get copies of the original bins.
    my @startBins = map { Bin->new_copy($_) } @$binList;
    # Compute the bins.
    my (undef, $bins) = Bin::Compute::ProcessScores(\@startBins, $score);
    # Analyze the bins.
    my $quality = Bin::Analyze::Quality($bins);
    my $scoreString = "[" . join(", ", @realScores) . "]";
    print "** Score computed for $scoreString is $quality. " . int(time - $start + 0.5) . " seconds for run.\n";
    my $stats = Bin::Analyze::Report($bins);
    print $oh "REPORT FOR $scoreString\n" . $stats->Show() . "\n\n";
    return $quality;
}

