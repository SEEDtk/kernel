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

=head1 Analyze Coding Kmers

    p3x-coding-kmers.pl [options]

This script inventories kmers in a set of genomes and classifies them according to whether they are non-coding, or coding
within a specific frame (+1,+2,+3,-1,-2,-3).  Only kmers that are wholly within or without a non-overlapping coding region
will be counted.  The intent is to identify kmers that can be used to identify coding regions in genomes that have not yet
been annotated.

=head2 Parameters

The positional parameter is the name of the output directory.  The output directory will contain
the database in JSON format in the file C<kmers.json>, and a summary of the useful kmers in
the file C<stats.tbl>.

The standard input should contain the IDs of the genomes to process.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  The column containing the input genome
IDs can be specified using L<P3Utils/col_options>.

The following additional command-line parameters are supported.

=over 4

=item verbose

Progress messages will be written to STDERR.

=item K

The DNA kmer size to use.  The default is C<15>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use KmerFramer;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('workDir', P3Utils::col_options(), P3Utils::ih_options(),
    ['verbose|debug|v', 'write progress messages to STDERR'],
    ['kmer|K=i', 'DNA kmer size', { default => 15 }],
    );
my $stats = Stats->new();
# Get the options.
my $K = $opt->kmer;
my $debug = $opt->verbose;
# Get the work directory name.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No output directory specified.";
} elsif (! -d $workDir) {
    print STDERR "Creating output directory $workDir.\n" if $debug;
    File::Copy::Recursive::pathmk($workDir) || die "Could not create working directory $workDir: $!";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Create the utility object.
my $kmerFramer = KmerFramer->new(stats => $stats, p3 => $p3, K => $K, debug => $debug);
# Loop through the input.
my ($batchCount, $genomeCount) = (0, 0);
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # First get the names of the genomes.  We rebuild the couplets to genome IDs only.
    $batchCount++;
    print STDERR "Retrieving genome batch $batchCount.\n" if $debug;
    my @couples = map { [$_->[0], [$_->[0]]] } @$couplets;
    my $nameResults = P3Utils::get_data_batch($p3, genome => [], ['genome_name'], \@couples);
    print STDERR scalar(@$nameResults) . " genomes found.\n" if $debug;
    # Loop through the genomes found.
    for my $genomeData (@$nameResults) {
        $stats->Add(genomesIn => 1);
        $genomeCount++;
        my ($genomeID, $genomeName) = @$genomeData;
        print STDERR "Processing genome $genomeCount: $genomeID $genomeName.\n" if $debug;
        # Get all the protein feature data for this genome.
        my $seqMap = $kmerFramer->SequenceMap($genomeID);
        # Now we run through the sequences, counting kmers.
        print STDERR "Retrieving sequences.\n" if $debug;
        my $seqList = P3Utils::get_data($p3, contig => [['eq', 'genome_id', $genomeID]],
                ['sequence_id', 'sequence']);
        while (my $seqData = pop @$seqList) {
            my ($seqID, $sequence) = @$seqData;
            print STDERR "Processing $seqID.\n" if $debug;
            $kmerFramer->CountKmers($sequence, $seqMap->{$seqID});
        }
    }
}
# Write the kmer database.
print STDERR "Creating output file in $workDir.\n" if $debug;
$kmerFramer->Store("$workDir/kmers.json");
# Compute the mean and standard deviation.
print STDERR "Computing statistical metrics.\n" if $debug;
my ($mean, $sdev, $kCount) = $kmerFramer->Metrics();
print STDERR "Mean is $mean with deviation $sdev over $kCount kmers.\n" if $debug;
# Create the distribution analysis.
my $dHash = $kmerFramer->Distribution($mean, $sdev);
# Write the brackets.
open(my $oh, '>', "$workDir/brackets.tbl") || die "Could not open brackets.tbl: $!";
P3Utils::print_cols(['bracket', 'count'], oh => $oh);
for my $bracket (sort { $a <=> $b } keys %$dHash) {
    P3Utils::print_cols([$bracket, $dHash->{$bracket}], oh => $oh);
}
close $oh; undef $oh;
print STDERR "All done.\n" . $stats->Show();

