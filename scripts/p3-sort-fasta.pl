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

=head1 Reorder a FASTA File Using a Template

    p3-sort-fasta.pl [options] genomeFile <fastaIn >fastaOut

This script reads a FASTA file and then outputs its records in the order indicated by a list of genomes.  The FASTA
file should have feature IDs for the sequence IDs.

=head2 Parameters

The positional parameter is the name of the genome list file.  The first genome ID in each line will be used.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use FastA;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('genomeFile', P3Utils::ih_options(),
    );
# Get the genome file.
my ($genomeFile) = @ARGV;
if (! $genomeFile) {
    die "No genome file specified."
} elsif (! -f $genomeFile) {
    die "Genome file $genomeFile not found or invalid.";
}
# Open the input file.
my $ih = P3Utils::ih($opt);
my $fastH = FastA->new($ih);
# We will stash our sequences in here, keyed by genome.
my %triples;
# Loop through the input.
while ($fastH->next()) {
    my $id = $fastH->id();
    if ($id =~ /(\d+\.\d+)/) {
        push @{$triples{$1}}, [$id, $fastH->comment(), $fastH->left()];
    }
}
# Now we have all the sequences.  Read the genome file and write the
# sequences in order.
open(my $gh, '<', $genomeFile) || die "Could not open genome file: $!";
while (! eof $gh) {
    my $line = <$gh>;
    if ($line =~ /(\d+\.\d+)/) {
        my $triplets = $triples{$1};
        for my $triplet (@$triplets) {
            print ">$triplet->[0] $triplet->[1]\n$triplet->[2]\n";
        }
    }
}