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
use Shrub;
use ScriptUtils;

=head1 Synthesize Metagenomes

    SynthesizeMeta.pl [ options ] genome1 genome2 ...

Synthesize metagenomes from a set of input genomes. This script uses L<art_illumina|http://www.niehs.nih.gov/research/resources/software/biostatistics/art/>
to generate the synthetic files. Currently, that product is not installed as part of the SEEDtk environment, so this
will need to be changed. 

=head2 Parameters

The positional parameters are a list of genome IDs.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item cov

Fold coverage used in simulating reads. The default is C<20>.

=item err

Fraction of standard Illumina 2000 error rate. The default is C<1.0>.

=item len

Simulated read length. The default is C<100>.

=item mean

Mean of paired end read fragment size. The default is C<300>.

=item stdev

Standard deviation of read fragment size. The default is C<10>.

=item dir

Name of the working directory. The output files (C<meta1.fq> and C<meta2.fq>) will be placed in here, as well as the
input file created from the input genome contigs.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('genome1 genome2 ... genomeN', 
        Shrub::script_options(),
        ['cov|c=f',     'fold coverage for simulated reads', { default => 20 }],
        ['err|e=f',     'error rate fraction', { default => 1.0 }],
        ['len|l=i',     'simulated read length', { default => 100 }],
        ['mean|m=f',    'mean of paired end read fragment size', { default => 300}],
        ['stdev|s=f',   'standard deviation of read fragment size', { default => 10 }],
        ['dir|d=s',     'working directory', { required => 1 }]
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
my @genomes = @ARGV;
# Check the working directory.
my $dir = $opt->dir;
die "Invalid working directory $dir." if (! -d $dir);
# Open the input file.
open(my $oh, '>', "$dir/input.fa") || die "Could not create input file: $!";
# Loop through the genomes.
for my $genome (@genomes) {
    my $dnaFile = $shrub->genome_fasta($genome);
    open(my $ih, '<', $dnaFile) || die "Could not open genome DNA file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        print $oh $line;
    }
}
# Close off the output file.
close $oh;
# Set up the parameters for the ART Illumina run.
my $qs = -10 * log($opt->err)/log(10);
my $art = "/home/fangfang/bin/art_illumina";
my $out = "$dir/input.fa";
my $len = $opt->len;
my $cov = $opt->cov;
my $mean = $opt->mean;
my $stdev = $opt->stdev;
my $cmd = "$art -ss HS20 -i $dir/input.fa -l $len -f $cov -p -m $mean -s $stdev -qs $qs -qs2 $qs -na -o $out &>art.log";
my $rc = system($cmd);
if ($rc) {
    die "FAILED $art: rc = $rc";
}
