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
use Shrub;
use Shrub::GTO;


=head1 Produce a Quality Report on Shrub Genomes

    shrub_report.pl [ options ] skreport cmreport

This script examines the information in two reports on the Shrub genomesto produce a quality report.

The columns in the report are as follows.

=over 4

=item 1

Genome source (currently CORE or PATRIC).

=item 2

Genome ID

=item 3

Genome name

=item 4

Number of contigs

=item 5

Number of DNA base pairs

=item 6

C<1> if the genome is mostly complete, else C<0>

=item 7

The contig N50 score.

=item 8

The contig N70 score.

=item 9

The contig N90 score.

=item 10

The genome ID (for compatibility with other reports).

=item 11

The genome name (for compatibility with other reports).

=item 12

SciKit coarse evaluation score (percent)

=item 13

SciKit fine evaluation score (percent)

=item 14

CheckM evaluation completeness (percent)

=item 15

CheckM evaluation contamination (percent)

=item 16

CheckM taxonomy classification

=back

The two input files work as follows.

=item cmreport

CheckM results file. The fourth and subsequent lines of this file contain space-delimited fields, the first of which is the genome
ID. The second field is the CheckM taxonomy classification. The thirteenth is the checkM completeness and the fourteenth is the
checkM contamination score.

=item skreport

SciKit evaluator results file, tab-delimited. The first column contains the genome ID, the third contains the coarse score
and the fourth contains the fine score.

=back

=head2 Parameters

The positional parameters are the name of the scikit report file and the checkm report file.

The command-line options are those in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('skreport cmreport', Shrub::script_options(),
        );
# Get the two report file names.
my ($skreport, $cmreport) = @ARGV;
if (! $skreport) {
    die "No SciKit report specified.";
} elsif (! -s $skreport) {
    die "SciKit report $skreport invalid or missing.";
} elsif (! $cmreport) {
    die "No CheckM report specified.";
} elsif (! -s $cmreport) {
    die "CheckM report $cmreport invalid or missing.";
}
# These hashes will contain the quality data for the reports.
my %quality;
# Read in the SciKit data.
print STDERR "Reading SciKit report.\n";
open(my $ih, '<', $skreport) || die "Could not open $skreport: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($genome, undef, $coarse, $fine) = split /\t/, $line;
    $quality{$genome} = [$coarse, $fine];
}
close $ih; undef $ih;
# Read in the Checkm data.
print STDERR "Reading CheckM report.\n";
open($ih, '<', $cmreport) || die "Could not open $cmreport: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    if ($line =~ /^\s+(\d+\.\d+)\s+(\S+).+?(\d+\.\d+)\s+(\d+\.\d+)\s+\d/) {
        my ($genome, $taxon, $score, $contam) = ($1, $2, $3, $4);
        my $list = $quality{$genome};
        if (! $list) {
            print STDERR "$genome in CheckM report but not SciKit.\n";
        } else {
            $quality{$genome} = [$list->[0], $list->[1], $score, $contam, $taxon];
        }
    }
}
# Discard incomplete genomes.
print STDERR "Verifying reports.\n";
for my $genome (keys %quality) {
    my $list = $quality{$genome};
    if (scalar @$list < 5) {
        print STDERR "$genome in SciKit report but not CheckM.\n";
        delete $quality{$genome};
    }
}
# Get the database.
print STDERR "Connecting to the database.\n";
my $shrub = Shrub->new_for_script($opt);
# Loop through the genomes, producing output.
my @genomes = $shrub->GetAll('Genome', '', [], 'id name core dna-size contigs');
print STDERR "Processing genomes.\n";
for my $genomeDesc (@genomes) {
    my ($genome, $name, $core, $dnaSize, $contigs) = @$genomeDesc;
    # Only proceed if we have quality data on this genome.
    my $list = $quality{$genome};
    if ($list) {
        # Get the quality data.
        my ($coarse, $fine, $checkm, $contam, $taxon) = @$list;
        # Compute the source.
        my $source = ($core ? 'CORE' : 'PATRIC');
        # Get the key metrics.
        my $gto = Shrub::GTO->new($shrub, $genome);
        my $metricsH = $gto->metrics();
        print join("\t", $source, $genome, $name, $contigs, $dnaSize,
                $metricsH->{complete}, $metricsH->{N50}, $metricsH->{N70},
                $metricsH->{N90}, $genome, $name, $coarse, $fine, $checkm, $contam,
                $taxon) . "\n";
    }
}
