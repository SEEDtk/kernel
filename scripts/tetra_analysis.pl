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
use TetraMap;
use TetraProfile;

=head1 Analyze Tetramer Vectors for Well-Behaved Genomes

    tetra_analysis.pl [ options ] 

This script runs through all the well-behaved genomes in the Shrub database and computes their tetramer profiles. For each genome, it will
compute the global tetramer profile and then the distance of each chunk from that profile. It will then output the mean, standard deviation,
minimum, maximum, and expected maximum (tolerance) for the distances. The distances are estimated from the dot-product of the normalized
profile vectors. 

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item chunkSize

The size of each DNA chunk to use when computing tetramers. The default is C<1000>.

=item input

If specified, a list of genome IDs to process. This list will be used instead of the list of well-behaved genomes.

=item dist

If specified, the distances will be true vector distance instead of 1 minus the dot product.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        ['input|i=s', 'name of file containing genome IDs to process'],
        ['chunkSize|l=i', 'chunk size for breaking up long contigs', { default => 1000 }],
        ['dist', 'use true distance']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Compute the distance type.
my $type = ($opt->dot ? 'dot' : 'dist');
print STDERR "Distance type is '$type'.\n";
# Get the list of genomes. For each, we need its name and DNA file.
my $dnaDir = $FIG_Config::shrub_dna;
my %genomes;
if ($opt->input) {
    # Here we have an input file.
    my $file = $opt->input;
    print STDERR "Reading genomes from $file.\n";
    open(my $ih, '<', $file) || die "Could not open input file $file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /(\d+\.\d+)/) {
            my $genome = $1;
            my ($gData) = $shrub->GetAll('Genome', 'Genome(id) = ?', [$genome], 'id name contig-file');
            my ($id, $name, $contigFile) = @$gData;
            $genomes{$id} = [$name, $contigFile];
        }
    }
} else {
    print STDERR "Reading genomes from database.\n";
    %genomes = map { $_->[0] => [$_->[1], $_->[2]] } $shrub->GetAll('Genome', 'Genome(well-behaved) = ?', [1], 'id name contig-file');
}
my $n = scalar keys %genomes;
print STDERR "$n genomes to process.\n";
# Get the chunk size and build the tetramap generator.
my $chunkSize = $opt->chunksize;
my $mapper = TetraMap->new();
# Output the header.
print join("\t", qw(genome name min max mean sdev tol)) . "\n";
# Loop through the genomes.
my $i = 0;
for my $genome (sort keys %genomes) {
    my ($name, $file) = @{$genomes{$genome}};
    print STDERR "Processing $genome $i of $n: $name\n";
    my $profile = TetraProfile->new($mapper, $chunkSize);
    $profile->ProcessFasta("$dnaDir/$file");
    print join("\t", $genome, $name, $profile->stats($type)) . "\n";
}
print STDERR "All done.\n";


## TODO process the input to produce the output.