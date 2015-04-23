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
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use File::Copy::Recursive;

=head1 Community Analysis Pipeline Script

        community_pipeline [-l MinLn] [-p MinPsc] [-s MinSim] -d DataDir -c Sample -g CalibrationGenomes

This script analyzes a set of contigs and separates them into clusters that likely represent DNA from a single
organism. In order to do this, the contigs are first compared against a specified set of reference genomes.
For each contig, a I<distance vector> is formed that indicates the similarity of that contig to each reference
genome. The distance between two contigs can be represented by a dot product of the vectors. The higher the
dot product, the closer the contigs. The contigs are then partitioned. Initially, each contig is in its own
partition. Partitions whose contigs are very close are merged, with the proviso that too many occurrences of
particular functional roles is an indicator the set has gotten too big.

=head2 Parameters

There are no positional parameters. The following command-line options are supported.

=over 4

=item minlen

=item maxpsc

=item minsum

=item normalize

=item data

The name of a master directory in which all the working and output files will be stored. If the directory
does not exist, it will be created. To save time when this script is invoked with differing tuning parameters,
intermediate files created here will be reused if they already exist.

=item sample

Name of a FASTA file containing the community sample to be analyzed. If a copy already exists in
the directory given by the C<--data> parameter, this parameter will be ignored.

=item samplename

The name given to this environmental sample. This is also the name of the subdirectory under the master data directory
(C<--data>) in which all output files are stored.

=item minkhits

The minimum number of kmer hits in a reference genome for that genome to be considered of interest in analyzing the
community sample.

=item refsf

The name of a file containing a list of the representative genomes. The community will be analyzed to find sets of
contigs each of which resembles one of these representative genomes. The file is tab-delimited, with the genome ID
in the first column and the genome name in the second.

=back

=head2 Output

All of the intermediate and output files are placed in a subdirectory of the master directory
having the same name as the incoming community sample (C<--samplename>). Thus, the master directory can contain
the output for numerous samples, so long as they all have different names.

=over 4

=item sample.fa

A copy of the FASTA file containing the contigs for the community sample.

=item repk.json

A JSON file containing kmer data relating to the representative genomes. If this file already exists, then
it will be re-used and the C<--refsf> parameter will be ignored.

=item ref.counts

A tab-delimited file containing the number of hits in the sample to indicative kmers from the representative
genomes, sorted by most hits to least hits. Each record in the file represents a single reference genomes for
which there were at least 5 kmer hits found. Each record contains three columns-- (0) the number of hits,
(1) the genome ID, and (2) the genome name.

=item RefD

A directory containing the reference genomes of interest to the sample. A reference genome is only considered
of interest if it has at least the number of kmer hits indicated by the C<--minkhits> parameter. Each genome
will have a subdirectory having the genome ID as its name. The subdirectory will contain the following files.

=over 8

=item genome.gto

A JSON-format L<GenomeTypeObject> for the genome.

=item reference.contigs

A FASTA file for the genome's contigs.

=item blast.out

The BLAST output from comparing the genome's contigs to the contigs in the community sample.

=back



=back

##TODO output description

=cut

my $opt =
  ScriptUtils::Opts( '',
                        [ 'minlen|l=i', 'minimum length for blast match to count', { default => 500 } ],
                        [ 'maxpsc|p=f', 'maximum pscore for blast match to count', { default => 1.0e-100 } ],
                        [ 'minsim|s=f', 'minimum % identity for condensing', { default => 7000 } ],
                        [ 'normalize', 'use normalized distances'],
                        [ 'data|d=s', 'Data Directory for Community Pipeline', { required => 1 } ],
                        [ 'sample|c=s','community DNA sample in fasta format', { } ],
                        [ 'samplename|n=s','environmental Sample Name', { required => 1 } ],
                        [ 'minkhits|k=i','minimum number hits to be a reference', { default => 400 } ],
                        [ 'refsf|r=s','File of potential reference genomes', { } ]
    );

my $dataD      = $opt->data;
my $sample     = $opt->sample;
my $sample_id  = $opt->samplename;
my $refsF      = $opt->refsf;
my $min_hits   = $opt->minkhits;
my $max_psc    = $opt->maxpsc;
my $min_len    = $opt->minlen;
my $min_sim    = $opt->minsim;


if (! -d $dataD) { mkdir($dataD,0777) || die "could not make $dataD" }

if (! -d "$dataD/$sample_id")
{
    mkdir("$dataD/$sample_id") || die "could not make $dataD/$sample_id";
}
my $dataS = "$dataD/$sample_id";
if (! -s "$dataS/sample.fa")
{
    File::Copy::Recursive::fcopy($sample, "$dataS/sample.fa");
}
if (! -s "$dataS/repk.json")
{
    if (! -s $refsF)
    {
        die "you need to use the --refs parameter to specify representative genomes";
    }
    &SeedUtils::run("compute_close_data > $dataS/repk.json < $refsF");
}
if (! -s "$dataS/ref.counts") {
    &SeedUtils::run("get_closest_representative_genomes -d $dataS/repk.json < $dataS/sample.fa > $dataS/ref.counts");
}
if (! -d "$dataS/RefD")
{
    &SeedUtils::run("construct_reference_genomes -c $dataS/sample.fa -m $min_hits -r $dataS/RefD < $dataS/ref.counts");
}
my $norm = $opt->normalize ? '-n' : '';
&SeedUtils::run("initial_estimate -r $dataS/RefD -c $dataS/ref.counts -l $min_len -p $max_psc -s $min_sim $norm -v $dataS/saved.sim.vecs > $dataS/bins");
&SeedUtils::run("summarize_bins < $dataS/bins > $dataS/bins.summary");
