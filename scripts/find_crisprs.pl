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
use crispr qw();
use SeedUtils qw();
use Shrub::Contigs;
use File::Copy::Recursive;
use Stats;

=head1 Search for CRISPR Arrays in Genomes

    find_crisprs.pl [ options ] outDir

This script will take as input a list of genome IDs and output CRISPR information for each identified genome.
For each such genome, a file will be created with the name I<XXXXXX.X>C<.json>, where I<XXXXXX.X> is the genome ID,
containing a JSON string encapsulating all the information needed to construct CRISPR features (repeat, array, spacer)
for that genome.

The JSON file will contain a list of I<CRISPR tuples>. Each tuple consists of (0) a location, (1) a consensus string
for the repeat region, (2) a list of repeat descriptors, and (3) a list of spacer descriptors. Each descriptor
is a 2-tuple consisting of (0) a location, and (1) a DNA sequence. Each location will be a 4-tuple containing (0) a contig
ID, (1) a beginning point, (2) a direction, and (3) a length.

=head2 Parameters

The single positional parameter should be the name of the output directory. If it does not exist, it will be created.

The standard input should be a tab-delimited file with genome IDs in the last column.

The command-line options are those found in L<Shrub/script_options> (to specify the database) and
L<ScriptUtils/ih_options> (to specify the standard input) plus the following.

=over 4

=item col

The (1-based) column index of the genome IDs in the input. The default is C<0>, indicating the last column.

=item missing

If specified, a genome will be scanned only if its output file does not already exist in the output directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['col|c=i', 'index (1-based) of genome ID column', { default => 0 }],
        ['missing|m', 'if specified, only new genomes will be processed']
        );
# Get the options.
my $missing = $opt->missing;
# Create a statistics object.
my $stats = Stats->new();
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (-f $outDir) {
    die "Invalid output directory specified.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Error creating $outDir: $!";
}
# Connect to the database.
print "Connecting to the database.\n";
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
print "Processing input.\n";
my $ih = ScriptUtils::IH($opt->input);
# Compute the input column index.
my $col = $opt->col - 1;
# Loop through the input.
while (! eof $ih) {
    # Locate the genome.
    my $line = <$ih>;
    $stats->Add(lineIn => 1);
    $line =~ s/[\r\n]+$//;
    my @fields = split /\t/, $line;
    my $genomeID = $fields[$col];
    if (! $genomeID) {
        $stats->Add(blankLine => 1);
    } elsif (-s "$outDir/$genomeID.json" && $missing) {
        print "Genome $genomeID already processed.\n";
        $stats->Add(genomeAlreadyRun => 1);
    } else {
        my $contigs = Shrub::Contigs->new($shrub, $genomeID);
        if (! $contigs) {
            print "Genome $genomeID not found.\n";
            $stats->Add(genomeNotFound => 1);
        } else {
            # Get the genome's sequence data.
            my @tuples = $contigs->tuples();
            my $contigCount = scalar @tuples;
            $stats->Add(contigsFound => $contigCount);
            # Run it through the CRISPR finder.
            print "Processing $genomeID with $contigCount contigs.\n";
            my $crisprData = crispr::find_crisprs(\@tuples);
            my $found = scalar @$crisprData;
            print "   $found arrays found.\n";
            $stats->Add(arraysFound => 1);
            if ($found) {
                $stats->Add(genomeHasArrays => 1);
            } else {
                $stats->Add(genomeNoArrays => 1);
            }
            # Write the CRISPR arrays to the output file.
            SeedUtils::write_encoded_object($crisprData, "$outDir/$genomeID.json");
        }
    }
}
print "All done:\n" . $stats->Show();

## TODO process the input to produce the output.