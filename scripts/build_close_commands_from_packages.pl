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

=head1 Create Representative-Genome Server Input File for Packages

    build_close_commands_from_packages.pl [ options ] packageDir

Loop through the genome packages in a directory. For each, extract the seed protein (Phenylalanyl tRNA synthetase alpha chain) and create a L<rep_server.pl>
command stream to find the closest genome.

=head2 Parameters

The positional parameter is the name of the genome package directory.

The command-line options are the following.

=over 4

=item closeness

The minimum degree of closeness desired. The default is C<50>.

=item seedProt

The name of the seed protein. The default is C<Phenylalanyl-tRNA synthetase alpha chain).

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir',
        ['closeness|K=i', 'minimum degree of similarity, in kmers', { default => 50 }],
        ['seedProt|seedprot|seed|p=s', 'functional assignment of seed protein', { default => 'Phenylalanyl-tRNA synthetase alpha chain' }]
        );
# Get the genome packages directory.
my ($packageDir) = @ARGV;
if (! $packageDir) {
    die "No package directory specified.";
} elsif (! -d $packageDir) {
    die "Invalid or missing package directory $packageDir.";
}
my $gHash = GPUtils::get_all($packageDir);
# Get the options.
my $closeness = $opt->closeness;
my $seedProt = $opt->seedprot;
# Loop through the genomes found.
for my $genome (sort keys %$gHash) {
    # Get the GTO.
    my $gto = GPUtils::gto_of($gHash, $genome);
    # Get the genome name.
    my $name = $gto->{scientific_name};
    # Extract the seed protein.
    my $featureList = GPUtils::role_to_features($gto, $seedProt);
    my $prot = $featureList->[0]{protein_translation};
    if (! $prot) {
        die "Seed protein not found in $genome: $name.";
    } else {
        print "# $genome $name\n";
        print "close_rep_seq $closeness $prot\n";
    }
}
