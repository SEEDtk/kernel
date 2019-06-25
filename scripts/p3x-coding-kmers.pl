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

The positional parameter is the name of the output directory.  In normal mode, the output directory will contain
the database in JSON format in the file C<kmers.json>, and a summary of the useful kmers in
the file C<stats.tbl>.  If the output directory is not specified, the program will run in file mode.  In this
case, the kmers will be written to the standard output as a two-column table, each record containing (0) a kmer, and
(1) a frame ID.

The standard input should contain the IDs of the genomes to process.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  The column containing the input genome
IDs can be specified using L<P3Utils/col_options>.

The following additional command-line parameters are supported.

=over 4

=item K

The DNA kmer size to use.  The default is C<15>.

=item shrub

If specified, the genomes are taken from a L<Shrub> database instead of PATRIC.  If this is the case, the Shrub
database is specified using L<Shrub/script_options>.

=item resume

If specified, then this is a restart of a partial job.  The datbase will be reloaded from C<kmers.json> in the
working directory.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Shrub;
use Stats;
use KmerFramer;
use File::Copy::Recursive;
use KmerFrameFiles;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('workDir', P3Utils::col_options(), P3Utils::ih_options(),
    Shrub::script_options(),
    ['shrub', 'use the Shrub database'],
    ['kmer|K=i', 'DNA kmer size', { default => 15 }],
    ['resume', 'resume after failed run'],
    );
my $stats = Stats->new();
# Get the options.
my $K = $opt->kmer;
my $fileMode = 0;
# Get the work directory name.
my ($workDir) = @ARGV;
if (! $workDir) {
    $fileMode = 1;
    print STDERR "File mode specified.\n";
} elsif (! -d $workDir) {
    print STDERR "Creating output directory $workDir.\n";
    File::Copy::Recursive::pathmk($workDir) || die "Could not create working directory $workDir: $!";
}
# Get access to the database.
my $p3;
if ($opt->shrub) {
    print STDERR "Connecting to Shrub.\n";
    $p3 = Shrub->new_for_script($opt);
} else {
    print STDERR "Connecting to PATRIC.\n";
    $p3 = P3DataAPI->new();
}
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Create the utility object.  How we do this depends on the resume flag.
my %options = (stats => $stats, p3 => $p3, debug => \*STDERR);
if ($opt->resume) {
    if ($fileMode) {
        die "Cannot resume in file mode.";
    } else {
        print STDERR "Loading saved results for resume.\n";
        $options{saved} = "$workDir/kmers.json";
    }
} else  {
    print STDERR "Creating new database.\n";
    $options{K} = $K;
}
my $kmerFramer = ($fileMode ? KmerFramerFiles->new(%options) : KmerFramer->new(%options));
# Loop through the input.
my ($batchCount, $genomeCount) = (0, 0);
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # First get the names of the genomes.  We rebuild the couplets to genome IDs only.
    $batchCount++;
    my $nameResults;
    eval {
        print STDERR "Retrieving genome batch $batchCount.\n";
        $nameResults = $kmerFramer->GenomeNames([map { $_->[0] } @$couplets]);
        print STDERR scalar(@$nameResults) . " genomes found.\n";
    };
    if ($@) {
        Checkpoint($kmerFramer, $@, $workDir);
    }
    # Loop through the genomes found.
    for my $genomeData (@$nameResults) {
        $stats->Add(genomesIn => 1);
        $genomeCount++;
        my ($genomeID, $genomeName) = @$genomeData;
        if ($kmerFramer->gCheck($genomeID)) {
            print STDERR "Already processed genome $genomeCount: $genomeID $genomeName.\n";
            $stats->Add(genomeSkipped => 1);
        } else {
            print STDERR "Processing genome $genomeCount: $genomeID $genomeName.\n";
            # We will fill these variables from the database.  The actual filling is protected so we can recover.
            my ($seqMap, $seqList);
            eval {
                # Get all the protein feature data for this genome.
                $seqMap = $kmerFramer->SequenceMap($genomeID);
                # Now we run through the sequences, counting kmers.
                print STDERR "Retrieving sequences.\n";
                $seqList = $kmerFramer->SequenceList($genomeID);
            };
            if ($@) {
                Checkpoint($kmerFramer, $@, $workDir);
            }
            print STDERR "Processing sequences.\n";
            while (my $seqData = pop @$seqList) {
                my ($seqID, $sequence) = @$seqData;
                $kmerFramer->CountKmers($sequence, $seqMap->{$seqID});
            }
            # Denote this genome is done.
            $kmerFramer->Record($genomeID);
            $stats->Add(genomeProcessed => 1);
        }
    }
}
# In normal mode, we do our output here.
if (! $fileMode) {
    # Write the kmer database.
    print STDERR "Creating output file in $workDir.\n";
    $kmerFramer->Save("$workDir/kmers.json");
    # Compute the mean and standard deviation.
    print STDERR "Computing statistical metrics.\n";
    my ($mean, $sdev, $kCount) = $kmerFramer->Metrics();
    print STDERR "Mean is $mean with deviation $sdev over $kCount kmers.\n";
    # Create the distribution analysis.
    my $dHash = $kmerFramer->Distribution($mean, $sdev);
    # Write the brackets.
    open(my $oh, '>', "$workDir/brackets.tbl") || die "Could not open brackets.tbl: $!";
    P3Utils::print_cols(['bracket', 'count'], oh => $oh);
    for my $bracket (sort { $a <=> $b } keys %$dHash) {
        P3Utils::print_cols([$bracket, $dHash->{$bracket}], oh => $oh);
    }
    close $oh; undef $oh;
}
print STDERR "All done.\n" . $stats->Show();


sub Checkpoint {
    my ($kmerFramer, $savedError, $workDir) = @_;
    if ($workDir) {
        # Here we failed retrieving data.  Checkpoint our results so far and percolate the error.
        print STDERR "ERROR retrieving genome data.  Saving progress.\n";
        $kmerFramer->Save("$workDir/kmers.json");
    }
    die "Fatal error: $savedError";
}