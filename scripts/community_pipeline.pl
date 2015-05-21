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

=item blast

The type of BLAST results to use when computing similarities to contigs in the
reference genomes. Both protein and nucleotide BLASTs are computed. If this
parameter is C<n>, the nucleotide results will be used. If it is C<p>, the
protein results will be used. A value of C<n> will cause more conservative
results. The default is C<n>.

=item minlen

The minimun length of a blast match necessary to consider that we have a correspondence between a
reference genome subsequence and a community sample contig subsequence. The default is C<500>. A
larger value will cause more conservative results.

=item maxpsc

The maximum permissible p-score (Poisson distribution rating) for a blast match to be acceptable.
This score is an indication of the degree to which a match may be random change. The default is
C<1e-100>. A lower value will cause more conservative results.

=item scoring

Scoring method for determining how to group contigs. The base algorithm groups the contigs
together by an abstract I<similarity score>. The similarity score is computed by comparing
I<scoring vectors>. The scoring vectors are computed by one of the following methods, as
determined by the value of this parameter.

=over 8

=item vector

The scoring vector contains the maximum percent identity score against each reference genome.

=item normvec

Same as C<vector>, but the vector is normalized to a unit length.

=item signal

The scoring vector contains the signal strength to the reference genome. The signal
strength between a sample contig and a reference genome is the sum of the percent identity
times the match length for each BLAST match, an indication of how much of the contig
matches the genome.

=item signalavg

The scoring vector contains the signal strength to the reference genome. The signal
strength between a sample contig and a reference genome is the mean of the percent identity
times the match length for each BLAST match, an indication of how much of the contig
matches the genome.

=back

=item compare

This specifies the algorithm for computing a similarity score from two scoring vectors.

=over 8

=item best

The score is an indication of how close together the two best scores in the vectors are. If the
best scores are at different positions, the score is 0.

=item bin

The score is C<1> if the scores are in the same relative order, else C<0>.

=item dist

The score indicates how close together the scores are at each vector position.

=item dot

The score is a dot product of the two vectors.

=back

=item minsim

The minimum similarity score for two community sample contigs to be considered for membership in
the same bin.

=item data

The name of a master directory in which all the working and output files will be stored. If the directory
does not exist, it will be created. To save time when this script is invoked with differing tuning parameters,
intermediate files created here will be reused if they already exist.

=item sample

Name of a FASTA file containing the community sample to be analyzed. If a copy already exists in
the directory given by the C<--data> parameter, this parameter will be ignored. The contig IDs
must be in the standard format for metagenomic assembly:
C<NODE_>I<id>C<_length_>I<contig_length>C<_cov_>I<coverage>C<_ID_>I<id2>.

=item samplename

The name given to this environmental sample. This is also the name of the subdirectory under the master data directory
(C<--data>) in which all output files are stored.

=item minkhits

The minimum number of kmer hits in a reference genome for that genome to be considered of interest in analyzing the
community sample.

=item refsf

The name of a file containing a list of the representative genomes. The community will be analyzed to find sets of
contigs each of which resembles one of these representative genomes. The file is tab-delimited, with the genome ID
in the first column and the genome name in the second. The default is C<representative.genomes> in the SEEDtk global
data directory.

=item univLimit

Maximum number of overlapping universal proteins allowed in a contig partition when the sample contigs are being grouped
together.

=item parmFile

If specified, the name of a file containing alternate parameter values. The file is space-delimited: each record begins
with a parameter name (long form) and is followed by one or more parameter values. All combinations of the parameter
values will be tried. The permissible parameters are C<blast>, C<minlen>, C<maxpsc>, C<minsim>, C<covgRatio>, C<univLimit>,
and C<normalize>. The output files (C<bins>, C<bins.summary>, and C<parms>) will be given numeric suffixes to distinguish them.

=item minCovg

The minimum coverage amount for a contig in the community sample. Only contigs with a coverage value greater than or equal
to this parameter will be included in the output. The default is C<0>, which means all are included.

=item expected

If specified, the name of a file containing the expected results. The file should be tab-delimited with two columns, the
first column being a contig ID and the second being a bin identifier. The results will be scored according to how well
they match the bin assignments in this file.

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

=item blast.out.dna

The BLAST output from comparing the genome's contigs to the contigs in the community sample.

=item blast.out.protein

The BLAST output from comparing the genome's proteins to the contigs in the community sample.

=back

=item saved.vecs

A tab-delimited file of similarity scores. The first record contains the blast type (C<p> or C<n>) followed
by a space and the scoring type. If these do not match the current values, the file is recomputed. Each additional
record in the file contains (0) a similarity score, (1) the ID of a source contig in the input
sample, (2) the ID of a second contig in the input sample.

=item merge.log

A log file of messages describing the binning choices made.

=item bins

A file describing in detail the partitions (bins) into which the contigs have been separated by the analysis
process. Each bin is described by a section terminated by a double-slash line (C<//>). First, the contigs from
the community sample that were assigned to the bin are listed. Each contig listing consists of

=over 8

=item 1

The contig ID, on a line by itself.

=item 2

One or more reference genome lines, sorted by the best match first. Each line contains (0) the similarity percentage
score between the reference genome and the contig, (1) the reference genome ID, and (2) the reference genome name.

=item 3

A blank line.

=back

After the contig section is a listing of universal roles found in the bin. For each universal role, there is a single
tab-delimited line consisting of the number of times the role was found followed by the role text.

=item bins.summary

The final output file. This file is divided into sections by bin, each section terminated by a double-slash line (C<//>).
For each bin, the file contains a summary of the universal role population followed by a summary of the reference genome
indications and a contig size count. These sections are separated by a blank line followed by a line of 8 hyphens (C<-------->).
The summary of universal role populations contains one line per role, tab-delimited, with each line containing
(0) a count of the number of occurrences of the role and (1) the role description. The summary of the reference genome
indications contains one line per genome, tab-delimited, with each line containing (0) the number of times the
reference genome was the closest match to a contig, (1) the total length of the contigs matched by the reference genome,
and (2) the name of the reference genome. The contig size count consists of the text C<total size of contigs>, an
equal sign (C<=>), and the total number of base pairs in all the contigs belonging to the bin.

=item parms

A tab-delimited file listing the major tuning parameters and their values.

=back

=cut

my $opt =
  ScriptUtils::Opts( '',
                        [ 'blast|b=s','blast type (p or n)',{ default => 'p'} ],
                        [ 'minlen|l=i', 'minimum length for blast match to count', { default => 500 } ],
                        [ 'maxpsc|p=f', 'maximum pscore for blast match to count', { default => 1.0e-100 } ],
                        [ 'minsim|s=f', 'minimum % identity for condensing', { default => 0.2 } ],
                        [ 'maxExpect|e=f', 'maximum E-value for BLASTing', { default => 1e-50 } ],
                        [ 'data|d=s', 'Data Directory for Community Pipeline', { required => 1 } ],
                        [ 'sample|c=s','community DNA sample in fasta format', { } ],
                        [ 'samplename|n=s','environmental Sample Name', { required => 1 } ],
                        [ 'minkhits|k=i','minimum number hits to be a reference', { default => 400 } ],
                        [ 'refsf|r=s','File of potential reference genomes', { default => "$FIG_Config::global/representative.genomes" } ],
                        [ 'covgRatio|cr=f', 'maximum acceptable coverage ratio for condensing', { default => 1.2 }],
                        [ 'univLimit|ul=i', 'maximum number of duplicate universal proteins in a set', { default => 2 }],
                        [ 'normalize=i', 'use normalized distances', { default => 0}],
                        [ 'parmFile=s', 'parameter specification file'],
                        [ 'minCovg|C=f', 'minimum coverage amount for a community sample contig', { default => 0 }],
                        [ 'scoring=s', 'scoring method', { default => 'vector' }],
                        [ 'compare=s', 'comparison method', { default => 'dotproduct' }],
                        [ 'expected=s', 'name of a file containing expected results'],
    );

my $blast_type = $opt->blast;
$blast_type    = ($blast_type =~ /^[pP]/) ? 'p' : 'n';

my $dataD      = $opt->data;
my $sample     = $opt->sample;
my $sample_id  = $opt->samplename;
my $refsF      = $opt->refsf;
my $maxE       = $opt->maxexpect;
my $min_hits   = $opt->minkhits;
my $expected   = $opt->expected;

my %parms;
$parms{blast}     = $blast_type;
$parms{maxpsc}    = $opt->maxpsc;
$parms{minlen}    = $opt->minlen;
$parms{minsim}    = $opt->minsim;
$parms{covgRatio} = $opt->covgratio;
$parms{univLimit} = $opt->univlimit;
$parms{minCovg}   = $opt->mincovg;
$parms{scoring}   = $opt->scoring;
$parms{compare}   = $opt->compare;

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
&SeedUtils::run("construct_reference_genomes -c $dataS/sample.fa -m $min_hits -e $maxE -r $dataS/RefD < $dataS/ref.counts");

if (! $opt->parmfile) {
    Process(\%parms);
} else {
    my %parmQ;
    open(PARMS, "<", $opt->parmfile) || die "Could not open parameter file: $!";
    while (defined (my $line = <PARMS>)) {
        chomp $line;
        my ($parm, @values) = split ' ', $line;
        if (! exists $parms{$parm}) {
            die "Invalid parameter name: $parm.";
        }
        $parmQ{$parm} = \@values;
    }
    close PARMS;
    my $runN = 1;
    my @keyQ = sort keys %parmQ;
    my @qPos = map { 0 } @keyQ;
    my @qLen = map { scalar @{$parmQ{$_}} } @keyQ;
    for my $parm (@keyQ) {
        $parms{$parm} = $parmQ{$parm}[0];
    }
    # We treat our positions in the parameter queues as a number with a variable
    # radix, and increment it each time. When a digit position overflows, we reset it to 0 and
    # increment the next position. So, if our parm value counts are 2,3,2 our numbers would go
    # 0,0,0 - 0,0,1 - 0,1,0 - 0,1,1 - 0,2,0 - 0,2,1 - 1,0,0 - etc.
    my $found = 1;
    while ($found) {
        print STDERR "Parameter run $runN.\n";
        Process(\%parms, $runN);
        $runN++;
        my $idx = $#keyQ;
        $found = 0;
        while (! $found && $idx >= 0) {
            $qPos[$idx]++;
            if ($qPos[$idx] >= $qLen[$idx]) {
                $qPos[$idx] = 0;
            } else {
                $found = 1;
            }
            $parms{$keyQ[$idx]} = $parmQ{$keyQ[$idx]}[$qPos[$idx]];
            $idx--;
        }
    }

}

sub Process {
    my ($parms, $suffix) = @_;
    my $realSuffix = (defined $suffix ? ".$suffix" : '');
    open(PARMS, ">$dataS/parms$realSuffix") || die "Could not open parms file: $!";
    for my $parm (sort keys %$parms) {
        print PARMS "$parm = $parms->{$parm}\n";
    }
    close PARMS;
    my $cmd = "initial_estimate -b $parms->{blast} -r $dataS/RefD -c $dataS/ref.counts -l $parms->{minlen} -p $parms->{maxpsc} " .
                    "-s $parms->{minsim} -v $dataS/saved.vecs -u $dataS/contig.to.uni --cr $parms->{covgRatio} --ul $parms->{univLimit} " .
                    "--minCovg $parms->{minCovg} --scoring $parms->{scoring} --compare $parms->{compare} " .
                    "--logFile $dataS/mergelog$realSuffix > $dataS/bins$realSuffix";
    &SeedUtils::run($cmd);
    my $expectation = ($opt->expected ? ("--expected " . $opt->expected) : "");
    &SeedUtils::run("summarize_bins -c $dataS/contig.to.uni $expectation < $dataS/bins$realSuffix > $dataS/bins.summary$realSuffix");
}

