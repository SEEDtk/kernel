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
use SeedUtils;
use P3DataAPI;
use File::Copy::Recursive;
use Stats;

=head1 Compute PATRIC Genome Quality

    patric_quality.pl [ options ]

This script will take as input a list of PATRIC genome IDs and will examine each to determine its quallity using
the SciKit quality measure. The quality number will be appended to the input file to produce the output, in the
manner of a standard pipeline script. That is, the standard input and output are both tab-delimited files, and each
output line will be an input line with the quality number at the end.

A minimum quality threshold can be specified, and in addition it is possible to request that the L<GenomeTypeObject>
files for the accepted genomes be saved in an output directory.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<ScriptUtils/ih_options> (which allows overriding the standard input) plus
the following.

=over 4

=item col

The index (1-based) of the input column containing the genome IDs. The default is C<0>, indicating the last column.

=item min

The minimum acceptable quality. Only acceptable genomes are included in the output. The default is C<0>.

=item save

The name of a directory where the L<GenomeTypeObject> files from acceptable genomes should be saved. The default is to
not save the genomes.

=item temp

The name of a temporary working directory. The default is the SEEDtk temporary directory. For best performance, this
directory should be on the same volume as the C<--save> directory (if one is specified). The actual working files will
be stored in a subdirectory having a name generated from the process ID.

=item complete

If specified, a completeness check will be performed before a genomes is checked for consistency. Incomplete genomes will
not be checked and will not be saved.

=item missing

If specified, only genomes without GTO files in the save directory will be checked.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        ['min|m=f',  'minimum acceptable quality score', { default => 0 }],
        ['save|S=s', 'save directory for GTO files'],
        ['col|c=i',  'genome ID column in input (1-based)', { default => 0 }],
        ['temp|t=s', 'temporary working directory', { default => $FIG_Config::temp }],
        ['missing',  'only check genomes without downloaded GTOs'],
        ['complete', 'skip incomplete genomes (based on an N70 check)'],
        );
# Create a statistics object.
my $stats = Stats->new();
# Figure out if we have a save directory.
my $saveDir = $opt->save;
if ($saveDir && ! -d $saveDir) {
    die "Invalid save directory $saveDir.";
}
# Get the control options.
my $missing = $opt->missing;
my $complete = $opt->complete;
# Get the minimum quality threshold.
my $min = $opt->min;
# Get the temporary directory.
my $tempDir = $opt->temp . "/p3q_$$";
File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
eval {
    # Create the results directory.
    my $resultDir = "$tempDir/results";
    File::Copy::Recursive::pathmk($resultDir) || die "Could not create $resultDir: $!";
    # Connect to PATRIC.
    my $p3 = P3DataAPI->new();
    # Open the input file.
    my $ih = ScriptUtils::IH($opt->input);
    # Loop through the input genomes.
    my @couplets;
    while (@couplets = ScriptUtils::get_couplets($ih, $opt->col, 1)) {
        my ($genomeID, $inputLine) = @{$couplets[0]};
        $stats->Add(lineIn => 1);
        print STDERR "Processing $genomeID.\n";
        # Check to see if it has already been processed.
        if ($saveDir && $missing && -s "$saveDir/$genomeID.gto") {
            # It's already processed, so we can skip it.
            $stats->Add(genomeSkipped => 1);
        } else {
            # Get the genome from PATRIC.
            my $gto = $p3->gto_of($genomeID);
            if (! $gto) {
                # Here the genome was not in PATRIC.
                $stats->Add(genomeNotFound => 1);
            } elsif ($gto->{domain} ne 'Bacteria' && $gto->{domain} ne 'Archaea') {
                # Here it is not prokaryotic.
                $stats->Add(genomeNotProk => 1);
            } elsif ($complete && ! $gto->is_complete) {
                # Here it is incomplete.
                $stats->Add(genomeIncomplete => 1);
            } else {
                # We have the genome and it's good, so save it to disk.
                my $tempGto = "$tempDir/$genomeID.gto";
                $gto->destroy_to_file($tempGto);
                # Clear the result directory.
                File::Copy::Recursive::pathrmdir($resultDir) || die "Could not clear $resultDir: $!";
                # Compute the quality.
                my $cmd = "gto_consistency $tempGto $resultDir $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $FIG_Config::global/roles.to.use";
                SeedUtils::run($cmd);
                # Read in the results.
                open(my $qh, "<$resultDir/evaluate.log") || die "Could not open $genomeID quality log: $!";
                my $quality;
                while (! eof $qh && ! defined $quality) {
                    my $line = <$qh>;
                    if ($line =~ /Consistency=\s+(\d+(?:\.\d+)?)%/) {
                        $quality = $1;
                    }
                }
                # Check for acceptability.
                if ($quality < $min) {
                    $stats->Add(genomeLowQuality => 1);
                } else {
                    print join("\t", @$inputLine, $quality) . "\n";
                    $stats->Add(genomeAccepted => 1);
                    if ($saveDir) {
                        File::Copy::Recursive::fmove($tempGto, "$saveDir/$genomeID.gto") ||
                            die "Could not save GTO file for $genomeID: $!";
                        $stats->Add(genomeSaved => 1);
                    }
                }
            }
        }
    }
};
if ($@) {
    print STDERR "FATAL ERROR: $@\n";
}
print STDERR "Cleaning up $tempDir.\n";
File::Copy::Recursive::pathrmdir($tempDir);
print STDERR "Statistics:\n" . $stats->Show();