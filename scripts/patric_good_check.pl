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
use P3DataAPI;
use GenomeTypeObject;
use GPUtils;
use Stats;

=head1 Look for Good PATRIC Genomes

    patric_good_check.pl [ options ]

This script locates good PATRIC genomes. A genome is considered good if it has good quality scores,
has exactly one occurrence of the seed protein, and the length of that occurrence is within normal limits.

The input file will be a list of genome IDs with an indication of which ones have good quality scores.

The standard output will be a list of genome IDs and names for the good genomes.

Frequent messages will be written to the standard error file, so it should be separate.

=head2 Parameters

The standard input should be a genome quality file. Such files are tab-delimited, and have the genome ID in the
second column and the genome name in the third. The 13th column contains the SciKit fine score and the 14th
column the CheckM score. The genome is of good quality if the SciKit fine score is 85 or more and the CheckM
score is 80 or more.

The command-line options are those found in L<ScriptUtils/ih_options> (for specifying the standard input) plus
the following.

=over 4

=item GTO

The name of a directory containing PATRIC GTO files. For performance reasons, this directory will be checked before
asking the PATRIC interface for the GTO. The GTO file name should be the genome ID followed by C<.gto>.

=item resume

If specified, must be the ID of a genome. It is presumed that an interrupted run is being resumed, and the genome ID
should be the first genome in the input file to process.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        ['GTO=s', 'GenomeTypeObject JSON file directory'],
        ['resume=s', 'resume interrupted run'],
        );
my $stats = Stats->new();
# This hash will map genome IDs to GTO file names.
my %gtos;
if ($opt->gto) {
    print STDERR "Scanning GTO directory.\n";
    my $gtoDir = $opt->gto;
    if (! -d $gtoDir) {
        die "Invalid GTO directory $gtoDir speciied.";
    } else {
        opendir(my $dh, $gtoDir) || die "Could not open GTO directory: $!";
        my @files = grep { $_ =~ /\.gto$/ } readdir $dh;
        closedir $dh;
        for my $file (@files) {
            $stats->Add(gtoDirFile => 1);
            if ($file =~ /^(\d+\.\d+)\.gto$/) {
                $gtos{$1} = "$gtoDir/$file";
                $stats->Add(gtoDirGenome => 1);
            }
        }
        print STDERR scalar(keys %gtos) . " genomes found in $gtoDir.\n";
    }
}
# Connect to PATRIC.
print STDERR "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# Check for a resume.
my $resumeGenome = $opt->resume;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
print STDERR "Reading input.\n";
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    $stats->Add(lineIn => 1);
    my @cols = split /\t/, $line;
    my $genome = $cols[1];
    if ($genome =~ /^\d+\.\d+$/) {
        # Check for a resume.
        if ($resumeGenome && $genome ne $resumeGenome) {
            $stats->Add(resumeSkip => 1);
        } else {
            $resumeGenome = '';
            my $count = $stats->Add(genomeIn => 1);
            if ($count % 100 == 0) {
                print STDERR "$count genomes read.\n";
            }
            my $name = $cols[2];
            if ($cols[12] >= 85 && $cols[13] >= 80) {
                $stats->Add(genomeGood => 1);
                # Get the GTO for this genome.
                my $gto;
                if ($gtos{$genome}) {
                    $gto = GenomeTypeObject->create_from_file($gtos{$genome});
                } else {
                    $gto = $p3->gto_of($genome);
                }
                if (! $gto) {
                    print STDERR "Genome $genome not found.\n";
                    $stats->Add(genomeNotFound => 1);
                } else {
                    my $isGood = GPUtils::good_seed($gto);
                    if ($isGood) {
                        $stats->Add(genomeKept => 1);
                        print "$genome\t$name\n";
                    }
                }
            }
        }
    }
}
print STDERR "All done.\n" . $stats->Show();
