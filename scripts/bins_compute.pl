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
use Bin::Compute;
use Bin::Score;
use Stats;

=head1 Build Bins From Community Contigs

    bins_compute [ options ] workDirectory <contigFile

This script reads contig information from the standard input and partitions them into bins based
on criteria relating to closest reference genomes, universal roles, coverage, and tetranucleotides.
The contigs are represented by L<Bin> objects.

=head2 Parameters

The single positional parameter is the name of a working directory to contain temporary and output files.

The command-line options are those found in L<ScriptUtils/ih_options> and L<Bin::Score/script_options>.

=head2 Input File

The input file contains one or more L<Bin> objects sequentially in the format described by L<Bin/Bin Exchange Format>.
Each represents a single contig from the input community.

=head3 Output Files

The output files are as follows, all in the working directory.

=over 4

=item bins.json

A file of L<Bin> objects in JSON format, suitable for reading by L<Bin/ReadBins>. These represent the computed bins.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('parms', ScriptUtils::ih_options(),
        Bin::Score::script_options(),
        );
# Create the scoring object.
my $score = Bin::Score->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the list of contigs. These are read in as bins.
print "Reading contigs from input.\n";
my $binList = Bin::ReadContigs($ih);
# Verify the working directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.\n";
}
# Compute the bins.
my ($stats, $bins) = Bin::Compute::Process($binList, $score);
# Write the resulting bins.
print "Writing bins.\n";
open(my $oh, ">", "$workDir/bins.json") || die "Could not open bins output file: $!";
for my $bin (@$bins) {
    $bin->Write($oh);
}
close $oh;
# Output the statistics.
print "All done.\n" . $stats->Show();
