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

=head1 Stitch Together Fasta and Genome Files in a RepServer

    fasta_stitch.pl [ options ] repServerDir

This script will put the genome names from the C<complete.genomes> file in the comment section of the associated C<6.1.1.20.fasta> file in a
L<rep_server.pl> directory.

=head2 Parameters

The positional parameter must be the name of the rep-server directory. The stitched FASTA file will be sent to the standard output.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('repServerDir');
# Get the input directory.
my ($repServerDir) = @ARGV;
if (! $repServerDir) {
    die "No input directory specified.";
} elsif (! -d $repServerDir) {
    die "Input directory $repServerDir not found or invalid.";
}
# Read in the genomes.
open(my $gh, "<$repServerDir/complete.genomes") || die "Could not open complete.genomes: $!";
my %genomes;
while (! eof $gh) {
    my $line = <$gh>;
    if ($line =~ /(\d+\.\d+)\s+(.+)/) {
        $genomes{$1} = $2;
    }
}
close $gh;
my $count = 0;
open(my $fh, "<$repServerDir/6.1.1.20.fasta") || die "Could not open 6.1.1.20.fasta: $!";
while (! eof $fh) {
    my $line = <$fh>;
    if ($line =~ /^>(.*?)(\d+\.\d+)(.+)?/) {
        my ($prefix, $genome, $comment) = ($1, $2, $3);
        $comment //= '';
        print ">$prefix$genome$comment $genomes{$genome}\n";
    } else {
        print $line;
    }
}
