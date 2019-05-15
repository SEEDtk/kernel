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

=head1 Scan and Fix a File of Representative-Genome Triples

    scan_triples.pl [options]

This script reads a set of representative-genome prototype files into memory, then it removes any for which there is more than
one instance per genome or the length is outside the mean range for the domain.  The mean range for the domain has to
be computed by this application.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  It should contain the output of
the following pipeline, where the desired role name is put in place of the string C<roleName>.

    p3-echo "roleName" | p3-role-features --attr=genome_id,aa_sequence | p3-join --key1=genome_id --only=taxon_lineage_ids,Score - patric.good.tbl

The input columns will be C<feature.genome_id>, C<feature.aa_sequence>, C<genome.taxon_lineage_ids>, and C<Score>.  We
need to reduce this to just the genome ID, score, and amino acid sequence, but only for genomes where the role is a valid
universal.

The command-line options are as follows.

=over 4

=item verbose

Display progress messages on STDERR.

=cut

use strict;
use P3Utils;
use Statistics::Descriptive;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('',P3Utils::ih_options(),
        ['verbose|debug|v', 'display progress messages on STDERR'],
    );
my $stats = Stats->new();
my $debug = $opt->verbose;
# Create the accumulator for the protein lengths.
my $prots = Statistics::Descriptive::Sparse->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my (undef, $cols) = P3Utils::find_headers($ih, roleProteins => 'genome_id', 'taxon_lineage_ids', 'Score', 'aa_sequence');
# Loop through the input, accumulating data.  For each genome, we need to know its protein string, its score, its domain,
# and whether or not it has too many proteins.
my %protHash;
my %scoreHash;
my %domainHash;
my %badHash;
print STDERR "Reading input file.\n" if $debug;
while (! eof $ih) {
    my ($genome, $lineage, $score, $seq) = P3Utils::get_cols($ih, $cols);
    $stats->Add(lineIn => 1);
    if ($scoreHash{$genome}) {
        # Here the genome has multiple proteins, so it is bad.
        $badHash{$genome} = 1;
        $stats->Add(duplicateProtein => 1);
    } else {
        # Store this genome's data.
        $protHash{$genome} = $seq;
        $scoreHash{$genome} = $score;
        if ($lineage =~ /::2157::/) {
            $domainHash{$genome} = 'A';
            $stats->Add(genomeA => 1);
        } elsif ($lineage =~ /::2::/) {
            $domainHash{$genome} = 'B';
            $stats->Add(genomeB => 1);
        } else {
            $badHash{$genome} = 1;
            $stats->Add(notProk => 1);
        }
    }
}
# Now we compute the mean and standard deviation.
print STDERR "Computing length limits.\n" if $debug;
my %statsH = ('A' => Statistics::Descriptive::Sparse->new(), 'B' => Statistics::Descriptive::Sparse->new());
for my $genome (keys %protHash) {
    if (! $badHash{$genome}) {
        my $domain = $domainHash{$genome};
        $statsH{$domain}->add_data(length($protHash{$genome}));
    }
}
# We now replace the statistics object with the tuple [minlength, maxlength];
for my $domain (keys %statsH) {
    my $statObject = $statsH{$domain};
    my $mean = $statObject->mean();
    my $sdev = $statObject->standard_deviation();
    my ($min, $max) = ($mean - 3*$sdev, $mean + 3*$sdev);
    print STDERR "Domain $domain has mean length $mean with range [$min, $max].\n" if $debug;
    $statsH{$domain} = [$min, $max];
}
# Now we output all the genomes with good lengths that are not bad.
print STDERR "Writing output.\n";
P3Utils::print_cols(['genome_id', 'score', 'aa_sequence']);
for my $genome (sort keys %protHash) {
    if ($badHash{$genome}) {
        $stats->Add(protTooMany => 1);
    } else {
        my $domain = $domainHash{$genome};
        my $seq = $protHash{$genome};
        my $len = length $seq;
        if ($len < $statsH{$domain}[0]) {
            $stats->Add(protTooShort => 1);
        } elsif ($len > $statsH{$domain}[1]) {
            $stats->Add(protTooLong => 1);
        } else {
            P3Utils::print_cols([$genome, $scoreHash{$genome}, $seq]);
        }
    }
}
print STDERR "All done.\n" . $stats->Show() if $debug;
