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
use Loader;

=head1 Create Bin FASTA Files

    bins_fasta.pl [ options ] workDir

This script takes the output of L<bins_generate.pl> and creates C<bin>I<X>C<.fa> files containing the contigs
for each bin. It is used when the bin FASTA files are needed outside of the standard L<bins_sample_pipeline.pl>
process (for example, when binning is done as a service inside another system).

=head2 Parameters

The single positional parameter is the name of the working directory for the binning sample. This directory must
contain the C<bins.json> file describing the bins and the C<sample.fasta> file containing the sample contigs.
These files are produced in the working directory specified on the L<bins_generate.pl> command.

If the FASTA files are to be submitted to BV-BRC for annotation, then an additional positional parameters is
required, specifying the workspace output path.  Each bin will have its own output directory inside the
folder specified.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir bv-brc-output-path');
# Get the loader object.
my $loader = Loader->new();
my $stats = $loader->stats;
# Compute the working directory.
my ($workDir, $outPath) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid or missing work directory $workDir.";
} elsif (! -s "$workDir/sample.fasta") {
    die "$workDir does not have a good sample.fasta file.";
} elsif (! -s "$workDir/bins.json") {
    die "$workDir does not have bins in it.";
}
# Compute the file names.
my $contigFastaFile = "$workDir/sample.fasta";
my $binJsonFile = "$workDir/bins.json";
# Read in the bins.
print "Reading bins from $binJsonFile.\n";
my $binList = Bin::ReadBins($binJsonFile);
# Loop through the bins, processing them one at a time. For each bin, we read the whole sample
# file to get the contigs. Only one bin's worth of contigs is kept in memory at a time, at the cost of
# a lot of extra reading.
my $binNum = 0;
for my $bin (@$binList) {
    my $taxonID = $bin->taxonID;
    my $name = $bin->name;
    $stats->Add(bins => 1);
    $binNum++;
    print "Processing bin $binNum - $taxonID: $name.\n";
    my %contigs = map { $_ => 1 } $bin->contigs;
    # Now we read the sample file and keep the contig triples.
    my @triples;
    my $ih = $loader->OpenFasta(sampleContigs => $contigFastaFile);
    my $fastaFile = "$workDir/bin$binNum.fa";
    open(my $oh, '>', $fastaFile) || die "Could not open FASTA output file for bin $binNum.";
    my $triple = $loader->GetLine(sampleContigs => $ih);
    while (defined $triple) {
        my $contigID = $triple->[0];
        if ($contigs{$contigID}) {
            $stats->Add(contigsKept => 1);
            push @triples, $triple;
            print $oh ">$triple->[0] $triple->[1]\n$triple->[2]\n";
        }
        $triple = $loader->GetLine(sampleContigs => $ih);
    }
    my $contigCount = scalar @triples;
    print "Found $contigCount contigs for bin $binNum.\n";
    if ($outPath) {
        # Here we want to submit the genome.
        my $rc = system("p3-submit-genome-annotation", -t => $taxonID, -n => $name, '--contigs-file' => $fastaFile,
                $outPath, "Bin.$binNum.$taxonID");
        print "Annotation request returned $rc.\n";
    }
}
print "All done.\n" . $stats->Show();

