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
use Contigs;
use BasicLocation;

=head1 Extract Locations From a FASTA File

    extract_markers.pl [ options ] contigs

This script creates a new FASTA file by removing specified locations from another FASTA file. The locations are
specified via the standard input.

=head2 Parameters

The single positional parameter is the name of the source FASTA file.

The command-line options are those found in L<ScriptUtils/ih_options> (to specify the input file with the locations) plus
the following.

=over 4

=item code

If specified, the output will be in the form of protein sequences using the genetic code specified for translation.

=back

The input file should be tab-delimited, with each line containing (0) the region label, (1) the region contig ID, (2) the
start location, (3) the strand, and (4) the length.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('contigs', ScriptUtils::ih_options(),
                ['code=i', 'genetic code for translation to protein (if omitted, DNA output)']
        );
# Get the parameter.
my ($contigs) = @ARGV;
if (! $contigs) {
    die "No contig FASTA file specified.";
} elsif (! -s $contigs) {
    die "Contig file $contigs missing or empty.";
}
# Get the genetic code, if any.
my $code = $opt->code;
# Read in the contigs.
my $contigObj = Contigs->new($contigs, genetic_code => $code);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the input, extracting the markers.
while (! eof $ih) {
    my ($label, $contig, $begin, $dir, $len) = ScriptUtils::get_line($ih);
    my $loc = BasicLocation->new($contig, $begin, $dir, $len);
    my $sequence;
    if ($code) {
        $sequence = $contigObj->xlate($loc);
    } else {
        $sequence = $contigObj->dna($loc);
    }
    print ">$label\n$sequence\n";
}
