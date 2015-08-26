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
use File::Copy::Recursive;
use Stats;

=head1 Process Single-Assembly Metagenomes

    bins_coverage [ options ] sampleFile outputDir

This script generates the coverage file (C<output.contigs2reads.txt>) for a single-sample metagenome. The coverage data is
taken from information in the FASTA file and used to create the output file.

=head2 Parameters

There are two positional parameters-- the name of the input FASTA file, and the name of the directory to contain the
prepared input files for the binning process.

The command-line options are as follows:

=over 4

=item keyword

Normally, the coverage is taken from the sequence ID in the FASTA file. Instead, the coverage can be coded as a keyword in the
sequence comments.  If so, the keyword name should be the value of this parameter. The coverage value must be connected to the
keyword name with an equal sign and the keyword must be only letters and digits.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleFile outputDir',
        ['keyword=s', 'keyword for coverage data in the sequence comments (if any)']
        );
# Get the positional parameters.
my ($sampleFile, $outputDir) = @ARGV;
if (! $sampleFile) {
    die "No input file name specified.";
} elsif (! -f $sampleFile) {
    die "Could not find input file $sampleFile.";
} elsif (! $outputDir) {
    die "No output directory specified.";
} elsif (! -d $outputDir) {
    die "Invalid output directory $outputDir.";
}
# Create the statistics object.
my $stats = Stats->new();
# Get the keyword parm.
my $keyword = $opt->keyword;
if ($keyword) {
    print "Using keyword mode.\n";
} else {
    print "Using standard contig ID mode.\n";
}
# Copy the input file to the output directory.
print "Copying contig file.\n";
File::Copy::Recursive::fcopy($sampleFile, "$outputDir/contigs.fasta");
# Open the input file.
open(my $ih, '<', $sampleFile) || die "Could not open contigs input file: $!";
# Open the output file.
open(my $oh, '>', "$outputDir/output.contigs2reads.txt") || die "Could not open output file: $!";
print "Computing coverage.\n";
# This will track bad coverage data.
my $errors = 0;
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(fastaLineIn => 1);
    # Is this an ID line?
    if ($line =~ /^>(\S+)\s+(.+)/) {
        # Yes. Get the ID and the comment.
        my ($contigID, $comment) = ($1, $2);
        $stats->Add(contigsFound => 1);
        # Look for the coverage.
        my $covg;
        if ($keyword) {
            if ($comment =~ /\b$keyword=([0-9.]+)/) {
                $covg = $1;
            } else {
                $stats->Add(missingKeyword => 1);
                $errors++;
            }
        } else {
            if ($contigID =~ /cov(?:erage|g)?_([0-9.]+)/) {
                $covg = $1;
            } else {
                $stats->Add(badContigID => 1);
                $errors++;
            }
        }
        if ($covg) {
            # We have coverage, produce the output line.
            print $oh "$contigID\t$covg\n";
            $stats->Add(contigOut => 1);
        }
    } else {
        # Not an ID line. Ignore it.
        $stats->Add(dataLineFound => 1);
    }
}
close $ih;
close $oh;
print "All done\n" . $stats->Show();
if ($errors) {
    print "\n** WARNING ** $errors errors found in input.\n";
}