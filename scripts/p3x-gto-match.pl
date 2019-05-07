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

=head1 Match GTOs Using Kmers

    p3x-gto-match.pl [options] gtoDir1 gtoDir2

This script will examine two directories of L<GenomeTypeObject> files, and attempt to match the ones that are most
similar.  A cross-reference of seed protein comparisons will be generated-- each GTO in the first directory will be
compared to each GTO in the second.  The bidirectional best hits based on protein kmers in common will be considered a
match.  The output file will list the IDs and names of the matching genomes.

=head2 Parameters

The positional parameters are the names of the directories containing the two sets of GTO files.

Additional command-line options are as follows.

=over 4

=item min

The minimum acceptable kmer match score for two genomes to be considered similar.  The default is C<100>.

=item verbose

Write progress messages to the standard error output.

=item protCheck

Paired genomes will be scored by the number of kmers in common between corresponding proteins.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use GeoGroup;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('gtoDir1 gtoDir2',
        ['min|m=i', 'minimum acceptable score for pairing', { default => 100 }],
        ['verbose|debug|v', 'write progress messages to the standard error output'],
        ['protCheck|prot', 'score genomes by protein BBHs']
);
my $stats = Stats->new();
# Get the options.
my $debug = $opt->verbose;
my $min = $opt->min;
my $details = $opt->protcheck;
# Get the input directories.
my ($gtoDir1, $gtoDir2) = @ARGV;
if (! $gtoDir1) {
    die "No input directories specified.";
} elsif (! -d $gtoDir1) {
    die "Input directory $gtoDir1 not found or invalid.";
} elsif (! $gtoDir2) {
    die "Only one input directory specified.";
} elsif (! -d $gtoDir2) {
    die "Input directory $gtoDir2 not found or invalid.";
}
# Compute the log file option.
my $logH = ($debug ? \*STDERR : undef);
print STDERR "Creating group 1 from $gtoDir1.\n" if $debug;
my $dLevel = ($details ? 2 : 0);
my $group1 = GeoGroup->new({ stats => $stats, logH => $logH, detail => $dLevel }, $gtoDir1);
print STDERR "Creating group 2 from $gtoDir2.\n" if $debug;
my $group2 = GeoGroup->new($group1->options, $gtoDir2);
print STDERR "Searching for bidirectional best hits.\n" if $debug;
# Compute the mapping between groups.
my ($pairs, $orphans1, $orphans2) = $group1->MapGroups($group2, $min);
my @headers = qw(genome1 name1 genome2 name2 seedScore);
push @headers, 'protScore' if $details;
P3Utils::print_cols(\@headers);
# Loop through the pairs, comparing the genomes.
for my $pair (@$pairs) {
    my ($g1, $g2) = @$pair;
    my $geo1 = $group1->geo($g1);
    my $geo2 = $group2->geo($g2);
    my @data = ($geo1->id, $geo1->name, $geo2->id, $geo2->name, $pair->[2]);
    if ($details) {
        print STDERR "Comparing $g1 to $g2.\n" if $debug;
        my $pHash1 = $geo1->protMap;
        my $pHash2 = $geo2->protMap;
        my @results = GEO::FindBBHs($pHash1, $pHash2, $min);
        my ($pCount, $o1Count, $o2Count) = map { scalar @$_ } @results;
        my $score = 0;
        if ($pCount) {
            $score = $pCount * 100 / ($pCount + ($o1Count + $o2Count) / 2);
        }
        push @data, $score;
    }
    P3Utils::print_cols(\@data);
}
# Write the orphan results.
my @trailer;
push @trailer, '' if $details;
print STDERR "Writing orphan output.\n" if $debug;
# Now do the orphans.
for my $genome (@$orphans1) {
    my $geo = $group1->geo($genome);
    P3Utils::print_cols([$geo->id, $geo->name, '', '', '', @trailer]);
    $stats->Add(orphan1 => 1);
}
for my $genome (@$orphans2) {
    my $geo = $group2->geo($genome);
    P3Utils::print_cols(['', '', $geo->id, $geo->name, '', @trailer]);
    $stats->Add(orphan2 => 1);
}
print STDERR "All done.\n" . $stats->Show() if $debug;
