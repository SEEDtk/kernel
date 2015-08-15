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

=head1 Search for Good Binning Parameters

    bins_genetic.pl [ options ] workDirectory

This program runs the community binning process multiple times to attempt to find the best parameter values,
using the L<AI::Genetic> genetic algorithm.

=head2 Parameters

The single positional parameter is the name of a working directory to contain temporary and output files.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item force

If specified, the scoring vectors will be recomputed even if the vector file already exists.

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

=item bins.json

A file of L<Bin> objects in JSON format, suitable for reading by L<Bin/ReadBins>. These represent the computed bins.

=item scores.tbl

A tab-delimited file. Each record contains two contig IDs followed by the elements of the scoring vector for the two
contigs. The vectors can be used to compute the final scores. If this file already exists, it will be reused.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDirectory', ScriptUtils::ih_options(),
        ['force', 'force recomputation of the scoring vectors'],
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
# Verify the working directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.\n";
}
# Compute the scoring vectors.
my $scoreVectors;
my $vectorFile = "$workDir/scores.tbl";
if ($opt->force || ! -f $vectorFile) {
    (undef, $scoreVectors) = Bin::Compute::ScoreVectors($binList);
    # Write out the scores.
    print "Saving score vectors to $vectorFile.\n";
    open(my $oh, ">", $vectorFile) || die "Could not open vector output file: $!";
    for my $scoreVector (@$scoreVectors) {
        my ($vector, $contig1, $contig2) = @$scoreVector;
        print $oh join("\t", $contig1, $contig2, @$vector) . "\n";
    }
} else {
    print "Reading score vectors from $vectorFile.\n";
    open(my $vh, "<", $vectorFile) || die "Could not open vector input file: $!";
    while (! eof $vh) {
        my $line = <$vh>;
        chomp $line;
        my ($contig1, $contig2, @vector) = split /\t/, $line;
        push @$scoreVectors, [\@vector, $contig1, $contig2];
    }
}
# Initialize the scoring parameter ranges.
$ga->init([[0, 0.5], [0, 0.5], [0, 0.5], [0, 0.5], [0, 1], [0.5, 1]]);
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
    # Create the scoring object.
    my $score = Bin::Score->new(@$scoreParms);
    # Compute the bins.
    my (undef, $bins) = Bin::Compute::ProcessScores($binList, $score, $scoreVectors);
    # Analyze the bins.
    my $quality = Bin::Analyze::Quality($bins);
    print "** Score computed for [" . join(", ", @$scoreParms) . "] is $quality. " . int(time - $start + 0.5) . " seconds for run.\n";
    return $quality;
}

