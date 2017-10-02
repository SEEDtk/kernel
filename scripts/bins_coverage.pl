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

=head1 Produce Metagenome Coverage Data

    bins_coverage [ options ] sampleFile outputDir

This script generates the coverage file (C<output.contigs2reads.txt>) for a metagenome. The coverage data is
taken from information in the FASTA file or assembly data files and used to create the output file.

=head2 Parameters

There are two positional parameters-- the name of the input FASTA file, and the name of the directory to contain the
prepared input files for the binning process.

The command-line options are as follows:

=over 4

=item standard

The coverage is encoded into the sequence ID from the FASTA file. Use this if the FASTA file was output from SPAdes.

=item keyword

The coverage coded as a keyword in the sequence comments.  The keyword name should be the value of this parameter.
The coverage value must be connected to the keyword name with an equal sign and the keyword must be only letters and digits.

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are retained. The default
is C<0>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are retained. The
default is C<0>.

=item covgFiles

If specified, a filename pattern for assembly data files. All the files matching this pattern will be parsed for
coverage data. The files are presumed to be tab-delimited, with the contig ID in the first column and the coverage
values in the second. The first line is treated as a header and discarded. The filename pattern must include the
directory path; otherwise the current directory is assumed.

=item noData

No coverage data is available. A dummy coverage of 50 is generated for each contig.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleFile outputDir',
        ['standard',     'parse coverage data from the sequence IDs'],
        ['keyword=s',    'keyword for coverage data in the sequence comments (if any)'],
        ['covgFiles=s',  'coverage data files (if any)'],
        ['noData',       'no coverage data available'],
        ['lenFilter=i',  'minimum contig length', { default => 0 }],
        ['covgFilter=f', 'minimum contig mean coverage', { default => 0}],
	['statistics-file=s', 'save statistics data to this file'],
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
# Figure out how we are computing coverage.
my $keyword = $opt->keyword;
my $covgFiles = $opt->covgfiles;
my $noData = $opt->nodata;
my $standard = $opt->standard;
if ($keyword) {
    if ($covgFiles) {
        die "covgFiles and keyword are mutually exclusive.";
    } elsif ($noData) {
        die "noData and keyword are mutually exclusive.";
    } elsif ($standard) {
        die "standard and keyword are mutually exclusive.";
    }
    print "Using keyword mode.\n";
} elsif ($covgFiles) {
    if ($noData) {
        die "covgFiles and noData are mutually exclusive.";
    } elsif ($standard) {
        die "covgFiles and standard are mutually exclusive.";
    }
    print "Using coverage file mode.\n";
} elsif ($standard) {
    if ($noData) {
        die "standard and noData are mutually exclusive.";
    }
    print "Coverage data will be parsed from sequence IDs.\n";
} elsif ($noData) {
    print "No coverage data available.\n";
} else {
    print "Coverage information will be inferred from sequence format.\n";
}
# This hash will contain coverage vectors if we are in file mode.
my %covgVectors;
# If in fact we are in file mode, we compute the coverage data here.
if ($covgFiles) {
    # Loop through the coverage files.
    my @files = grep { -f $_ } glob($covgFiles);
    my $pos = 0;
    for my $file (@files) {
        $stats->Add(covgFileIn => 1);
        open(my $ih, "<$file") || die "Could not open coverage file $file: $!";
        # Discard the header record.
        my $line = <$ih>;
        # Run through the other records, putting the coverage values at the current vector position.
        while (! eof $ih) {
            $line = <$ih>;
            chomp $line;
            $stats->Add(covgLineIn => 1);
            my ($contig, $covgValue) = split /\t/, $line;
            $covgVectors{$contig}[$pos] = $covgValue;
        }
        # Move to the next vector position.
        $pos++;
    }
    if (! $pos) {
        die "No coverage files found matching pattern $covgFiles.";
    }
    # Insure all the gaps are filled in.
    for my $contig (keys %covgVectors) {
        my $covgVector = $covgVectors{$contig};
        for (my $i = 0; $i < $pos; $i++) {
            if (! defined $covgVector->[$i]) {
                $stats->Add(covgGap => 1);
                $covgVector->[$i] = 0;
            }
        }
    }
}
# Open the input file.
open(my $ih, '<', $sampleFile) || die "Could not open contigs input file: $!";
# Open the coverage output file.
open(my $oh, '>', "$outputDir/output.contigs2reads.txt") || die "Could not open coverage output file: $!";
# Write the header.
print $oh "Contig ID\tCoverages\n";
# Open the FASTA output file.
open(my $fh, '>', "$outputDir/contigs.fasta") || die "Could not open FASTA output file: $!";
print "Reading FASTA input.\n";
# This will track bad coverage data.
my $errors = 0;
# These will track the coverage of the current contig.
my ($covg, $covgMean);
# These will track the ID, comment, and sequence of the current contig.
my ($contigID, $comment, $seq);
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(fastaLineIn => 1);
    # Is this an ID line?
    if ($line =~ /^>(\S+)(?:\s+(.+))?/) {
        # Yes. Get the ID and the comment.
        my ($newContigID, $newComment) = ($1, $2);
        # Insure we have a comment.
        $newComment //= '';
        $stats->Add(contigsFound => 1);
        # Process the previous contig.
        ProcessContig($fh, $stats, $contigID, $comment, $seq, $covgMean);
        # Initialize for the new contig.
        ($contigID, $comment, $seq) = ($newContigID, $newComment, "");
        # Look for the coverage.
        ($covg, $covgMean) = (0, 0);
        if ($covgFiles) {
            $covg = $covgVectors{$contigID};
            if (! defined $covg) {
                $stats->Add(contigNotFound => 1);
                $covg = 0;
                $errors++;
            } else {
                map { $covgMean += $_ } @$covg;
                $covg = join("\t", @$covg);
            }
        } else {
            if (! $keyword && ! $standard && ! $noData) {
                # Here we must compute the type of sequence label.
                if ($comment =~ /\b(covg|coverage|cov|multi)=[0-9.]+/) {
                    $keyword = $1;
                    print "Keyword \"$keyword\" selected.\n";
                } elsif ($contigID =~ /(?:coverage|covg|cov)_[0-9.]+/) {
                    $standard = 1;
                    print "Coverage will be parsed from sequence IDs.\n";
                } else {
                    $noData = 1;
                    print "No coverage data available.\n";
                }
            }
            if ($keyword) {
                if ($comment =~ /\b$keyword=([0-9.]+)/) {
                    $covg = $1;
                    $covgMean = $covg;
                } else {
                    $stats->Add(missingKeyword => 1);
                    $errors++;
                }
            } elsif ($noData) {
                $covg = 50;
                $covgMean = 50;
            } else {
                if ($contigID =~ /(?:coverage|covg|cov)_([0-9.]+)/) {
                    $covg = $1;
                    $covgMean = $covg;
                } else {
                    $stats->Add(badContigID => 1);
                    $errors++;
                }
            }
        }
        if ($covg) {
            # We have coverage, produce the coverage output line.
            print $oh "$contigID\t$covg\n";
            $stats->Add(coverageOut => 1);
            # Produce the statistical analysis of the coverage.
            my $covgCat;
            if ($covgMean >= 20) {
                $covgCat = '2X';
            } elsif ($covgMean >= 10) {
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
ProcessContig($fh, $stats, $contigID, $comment, $seq, $covgMean);
close $ih;
close $oh;
close $fh;
print "All done\n" . $stats->Show();
if ($opt->statistics_file)
{
    if (open(S, ">", $opt->statistics_file))
    {
	print S $stats->Show();
	close(S);
    }
    else
    {
	warn "Cannot write statistics file " . $opt->statistics_file . ": $!";
    }
}
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
