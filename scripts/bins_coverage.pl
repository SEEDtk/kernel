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

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are retained. The default
is C<0>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are retained. The
default is C<0>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleFile outputDir',
        ['keyword=s', 'keyword for coverage data in the sequence comments (if any)'],
                ['lenFilter=i',  'minimum contig length', { default => 0 }],
                ['covgFilter=f',  'minimum contig mean coverage', { default => 0}],
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
# Open the input file.
open(my $ih, '<', $sampleFile) || die "Could not open contigs input file: $!";
# Open the coverage output file.
open(my $oh, '>', "$outputDir/output.contigs2reads.txt") || die "Could not open coverage output file: $!";
# Open the FASTA output file.
open(my $fh, '>', "$outputDir/contigs.fasta") || die "Could not open FASTA output file: $!";
print "Computing coverage.\n";
# This will track bad coverage data.
my $errors = 0;
# These will track the coverage of the current contig.
my ($covg);
# These will track the ID, comment, and sequence of the current contig.
my ($contigID, $comment, $seq);
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(fastaLineIn => 1);
    # Is this an ID line?
    if ($line =~ /^>(\S+)\s+(.+)/) {
        # Yes. Get the ID and the comment.
        my ($newContigID, $newComment) = ($1, $2);
        $stats->Add(contigsFound => 1);
        # Process the previous contig.
        ProcessContig($fh, $stats, $contigID, $comment, $seq, $covg);
        # Initialize for the new contig.
        ($contigID, $comment, $seq) = ($newContigID, $newComment, "");
        # Look for the coverage.
        $covg = 0;
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
            # We have coverage, produce the coverage output line.
            print $oh "$contigID\t$covg\n";
            $stats->Add(coverageOut => 1);
            # Produce the statistical analysis of the coverage.
            my $covgCat;
            if ($covg >= 20) {
                $covgCat = '2X';
            } elsif ($covg >= 10) {
                $covgCat = '1X';
            } else {
                $covgCat = '0X';
            }
            $stats->Add("covg-$covgCat" => 1);
        }
    } else {
        # Not an ID line. Accumulate the sequence.
        chomp $line;
        $seq .= $line;
        $stats->Add(dataLineFound => 1);
    }
}
# Process the last contig.
ProcessContig($fh, $stats, $contigID, $comment, $seq, $covg);
close $ih;
close $oh;
close $fh;
print "All done\n" . $stats->Show();
if ($errors) {
    print "\n** WARNING ** $errors errors found in input.\n";
}


sub ProcessContig {
    my ($fh, $stats, $contigID, $comment, $seq, $covg) = @_;
    # Only process if we have a contig ID.
    if ($contigID) {
        my $len = length($seq);
        # Do we want this contig?
        if ($covg < $opt->covgfilter) {
            $stats->Add(contigBadCoverage => 1);
        } elsif ($len < $opt->lenfilter) {
            $stats->Add(contigTooShort => 1);
        } else {
            # Yes. Write the contig to the FASTA output.
            print $fh ">$contigID $comment\n$seq\n";
            $stats->Add(contigOut => 1);
        }
        # Record the length statistics.
        my $lenCat;
        if ($len < 1000) {
            $lenCat = '0000K';
        } else {
            my $lenThing = 1000000;
            my $zeroes = "";
            my $xes = "XXXK";
            while ($len < $lenThing) {
                $lenThing /= 10;
                $zeroes .= "0";
                $xes = substr($xes, 1);
            }
            my $cat = int($len / $lenThing);
            $lenCat = "$zeroes$cat$xes";
        }
        $stats->Add("contigLen-$lenCat" => 1);
        $stats->Add(letters => $len);
    }
}