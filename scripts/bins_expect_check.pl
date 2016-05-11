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

=head1 Analyze Bins Found Against Expectations

    bins_patric_check.pl [ options ] binDir

This script compares the species in the expected bins to the species actually found and the species currently
in the PhenTranSyntAlph database. There are three files of interest.

=over 4

=item 1

The C<ref.genomes.scores.tbl> file in the binning directory contains the reference genomes that matched
contigs in the binned metagenome. This is the I<found file>.

=item 2

The I<expectation file> specified on the command line contains the species expected in the metagenome.
This is the standard input file to the command.

=item 3

The C<PhenTranSyntAlph.fa> file in the Global data directory is a FASTA file of the primer proteins in the
PATRIC database. The comments in this file are the IDs and names of the corresponding genomes. This is the
I<primer file>.

=back

We focus on the first two words in each genome name, which are the genus and species. We are interested in
species in the found file that are not in the expectation file, species in the expectation file that are
not in the found file or the primer file, and species in the expectation file that are in the primer file
but not in the found file.

=head2 Parameters

The single positional parameter is the name of the binning directory.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item col

Input column (1-based) from the expectation file containing the genome names. The default is C<1>.

=item dcol

Input column (1-based) from the expectation file containing the abundance depth. The default is C<3>.

=item primer

The name of the primer file. The default is C<PhenTrnaSyntAlph.fa> in the global data directory.

=item genus

Perform the analysis on a genus rather than a genus/species basis.

=item mindepth

Minimum acceptable depth for an expected genome. The default is C<5>.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir', ScriptUtils::ih_options(),
        ['col|c=i', 'genome name input column', { default => 1 }],
        ['primer|p=s', 'primer file name', { default => "$FIG_Config::global/PhenTrnaSyntAlph.fa" }],
        ['genus|G', 'classify genomes by genus instead of genus and species'],
        ['dcol|d=i', 'abundance depth input column', { default => 3 }],
        ['mindepth|m=f', 'minimum depth for a genome to be considered expected', { default => 5 }],
        );
my $stats = Stats->new();
# Get the binning directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "Missing bin directory.";
} elsif (! -d $binDir || ! -s "$binDir/ref.genomes.scores.tbl") {
    die "Invalid bin directory $binDir.";
}
# Determine the number of words used for a genome category.
my $words = ($opt->genus ? 1 : 2);
# Read the primer file to get the primer groups.
my $primerFile = $opt->primer;
my %primersH;
print "Reading primer file $primerFile.\n";
open(my $ph, '<', $primerFile) || die "Could not open primer file: $!";
while (! eof $ph) {
    my $line = <$ph>;
    $stats->Add(primerLineIn => 1);
    if ($line =~ /^>\S+\s+\d+\.\d+\t(.+)/) {
        my $cat = name_to_cat($1);
        $stats->Add(primerHeaderIn => 1);
        if ($primersH{$cat}) {
            $primersH{$cat}++;
        } else {
            $primersH{$cat} = 1;
            $stats->Add(primerCatFound => 1);
        }
    }
}
close $ph; undef $ph;
# Now read the found file.
my $foundFile = "$binDir/ref.genomes.scores.tbl";
my %foundH;
print "Reading found file $foundFile.\n";
open($ph, '<', $foundFile) || die "Could not open found file: $!";
while (! eof $ph) {
    my $line = <$ph>;
    $stats->Add(foundLineIn => 1);
    if ($line =~ /\t([^\t]+)$/) {
        my $cat = name_to_cat($1);
        $stats->Add(foundCatLine => 1);
        if ($foundH{$cat}) {
            $foundH{$cat}++;
        } else {
            $foundH{$cat} = 1;
            $stats->Add(foundCatFound => 1);
        }
    }
}
close $ph; undef $ph;
# Open the input file.
print "Reading input for expected categories.\n";
my $ih = ScriptUtils::IH($opt->input);
# Discard the header line.
my $line = <$ih>;
# Now read it for the expected categories and depths.
my %expectedH;
my $col = $opt->col - 1;
my $dcol = $opt->dcol - 1;
my $mindepth = $opt->mindepth;
while (! eof $ih) {
    $line = <$ih>;
    $line =~ s/[\r\n]+$//; # chomp doesn't always work with stdin (ow!)
    $stats->Add(expectedLineIn => 1);
    my @cols = split /\t/, $line;
    my $name = $cols[$col];
    my $depth = $cols[$dcol];
    # Only proceed if we have good depth.
    if ($depth < $mindepth) {
        $stats->Add(expectedLowDepth => 1);
    } else {
        $stats->Add(expectedKept => 1);
        my $cat = name_to_cat($name);
        if ($expectedH{$cat}) {
            my $tuple = $expectedH{$cat};
            $tuple->[0]++;
            $tuple->[1] += $depth;
        } else {
            $expectedH{$cat} = [1, $depth];
            $stats->Add(expectedCatFound => 1);
        }
    }
}
# Convert the depth to an average.
for my $cat (keys %expectedH) {
    $expectedH{$cat}[1] /= $expectedH{$cat}[0];
}
# Loop through the expected categories.
print "Expected categories\n\nCat\tExpected\tDepth\tPrimer\tFound\n";
my @cats = sort { $expectedH{$b}[1] <=> $expectedH{$a}[1] } keys %expectedH; 
for my $cat (@cats) {
    my $primer = $primersH{$cat} // '';
    my $found = $foundH{$cat} // '';
    my ($count, $depth) = @{$expectedH{$cat}};
    print "$cat\t$count\t$depth\t$primer\t$found\n";
    if ($primer && ! $found) {
        $stats->Add('primers-not-found' => 1);
    }
}
print "\n\n";
# Loop through the found categories.
print "Unexpected categories\nCat\tFound\n";
for my $cat (sort keys %foundH) {
    my $expected = $expectedH{$cat};
    if (! $expected) {
        $stats->Add('found-not-expected' => 1);
        print "$cat\t$foundH{$cat}\n";
    }
}
print "\n\n";
# All done.
print "All done.\n" . $stats->Show();

# Convert a genome name to a category.
sub name_to_cat {
    my ($name) = @_;
    my @parts = split ' ', $name;
    splice @parts, $words;
    return join(' ', @parts);
}
