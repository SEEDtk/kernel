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
use GPUtils;
use Stats;

=head1 Extract SSU RNA sequences from Genome Packages

    package_rnas.pl [ options ] inDir

This script will run through all the genome packages in a directory and extract the DNA for the SSU rRNA features
into a FASTA file. The FASTA file will be produced on the standard output. Statistics and progress messages will
be sent to the standard error output, so this should be kept separate.

=head2 Parameters

The positional parameter is the name of the input directory.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir',
        );
# Create the statistics object.
my $stats = Stats->new();
# Validate the input directory.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Input directory $inDir missing or invalid.";
}
# Read the list of genomes.
my $genomeHash = GPUtils::get_all($inDir);
my $total = scalar(keys %$genomeHash);
print STDERR "$total genome packages found in input directory.\n";
# Loop through them.
my $count = 0;
for my $genome (sort keys %$genomeHash) {
    $stats->Add(genomeFound => 1);
    # Get the GTO.
    my $gto = GPUtils::gto_of($genomeHash, $genome);
    # Find the SSU role.
    my $featureList = GPUtils::role_to_features($gto, 'SSU rRNA');
    my $found = scalar @$featureList;
    print STDERR "$found SSU ribosomal RNAs found in $genome.\n";
    if (! $found) {
        $stats->Add(noRnaFound => 1);
    } else {
        # Loop through the features found.
        for my $feature (@$featureList) {
            my $fid = $feature->{id};
            my $dna = $gto->get_feature_dna($feature);
            print ">$fid\n$dna\n";
            $stats->Add(featureOut => 1);
        }
    }
}
print STDERR "All done.\n" . $stats->Show();
