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
use KmerDb;
use FastQ;

=head1 Check Samples for Likely Bins

    sample_kmer_check.pl [ options ] sampleDir kmerJson

This script examines the samples in a directory and marks them for likely bins.  Each read of the sample is examined for kmer
hits against a kmer database formed from seed protein DNA.  The maximum number of kmer hits for a single seed protein is
computed for the read, and it is considered likely if it has five or more. It is considered distant if it has less than ten.
The number of likely reads is output along with the number of likely and distant reads for each sample.

=head2 Parameters

The positional parameters are the name of the input sample directory and the name of the kmer database file. The kmer database
file should be the output file from L<seed_kmer_db.pl>.

The command-line options are the following.

=over 4

=item minHits

The minimum number of hits for a read to be considered likely. The default is C<5>.

=item maxHits

The maximum number of hits for a read to be considered distant. The default is C<10>.

=item verbose

Write status messages to the standard error output.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir kmerJson',
        ['minHits|min|m=i', 'minimum number of hits for likely reads', { default => 5 }],
        ['maxHits|max|M=i', 'maximum number of hits for distance reads', { default => 10 }],
        ['verbose|debug|v', 'write status messages to STDERR'],
        );
my $stats = Stats->new();
# Get the debug flag.
my $debug = $opt->verbose;
# Get the hit limits.
my $min = $opt->minhits;
my $max = $opt->maxhits;
# Get the file and directory names.
my ($sampleDir, $kmerJson) = @ARGV;
if (! $sampleDir) {
    die "No sample directory specified.";
} elsif (! -d $sampleDir) {
    die "Sample directory $sampleDir invalid or missing.";
} elsif (! $kmerJson) {
    die "No kmer file specified.";
} elsif (! -s $kmerJson) {
    die "Invalid or missing kmer file $kmerJson.";
}
# Load the kmer database.
print STDERR "Loading $kmerJson.\n" if $debug;
my $kmerdb = KmerDb->new(json => $kmerJson);
# Get the list of samples.
print STDERR "Reading $sampleDir.\n" if $debug;
opendir(my $dh, $sampleDir) || die "Could not open $sampleDir: $!";
my @samples = grep { -s ("$sampleDir/$_/$_" . '_1.fastq') } readdir $dh;
closedir $dh;
my $totSamples = scalar @samples;
print STDERR "$totSamples samples found in directory.\n" if $debug;
# Write the output headers.
print "Sample\tlikely\tl/d\tmax hits\tnulls\n";
# Process each sample.
for my $sample (@samples) {
    $stats->Add(sampleIn => 1);
    # Open the paired fastQ files.
    my @fq = map { "$sampleDir/$sample/${sample}_$_.fastq" } qw(1 2);
    my $fqh = FastQ->new(@fq);
    print STDERR "Processing $sample.\n" if $debug;
    # Loop through the reads, counting the likely ones and the likely-distant ones.
    my ($likely, $l_d, $maxHits, $nulls) = (0, 0, 0, 0);
    while ($fqh->next) {
        # Process the sequences.
        $stats->Add(readIn => 1);
        my %hits;
        for my $dna ($fqh->left, $fqh->right) {
            $kmerdb->count_hits($dna, \%hits);
        }
        # Get the maximum hit count.
        my ($best) = sort { $b <=> $a } values %hits;
        # Record the hit in the statistics.
        if (! $best) {
            $stats->Add(hitsFoundNone => 1);
            $best = 0;
            $nulls++;
        } elsif ($best < 10) {
            $stats->Add("hitsFound$best" => 1);
        } else {
            my $cat = int($best / 10) . "X";
            $stats->Add("hitsFound$cat" => 1);
        }
        if ($best >= $min) {
            $likely++;
            $stats->Add(likelyRead => 1);
            if ($best <= $max) {
                $l_d++;
                $stats->Add(likelyDistRead => 1);
            }
        }
        if ($best > $maxHits) {
            $maxHits = $best;
        }
    }
    # Write the counts.
    print join("\t", $sample, $likely, $l_d, $maxHits, $nulls) . "\n";
}
print STDERR "All done.\n" . $stats->Show() if $debug;