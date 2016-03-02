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
use Stats;
use SeedUtils;
use Loader;

=head1 Analyze Community Contigs to Determine Filtering Effect

    bins_filter_check.pl [ options ] sampleDir workDir

Process the crAss output for community samples and determine the number of bins that would survive a
filtering check.

=head2 Parameters

There are two positional parameters

=over 4

=item 1

The name of a directory containing the sample. The sample's contigs must be in a file called C<contigs.fasta> and
the vector file in a file called C<output.contigs2reads.txt> in this directory.

=item 2

The name of the directory to contain the output and intermediate files. If this parameter is omitted, the input
directory is assumed.

=back

The command-line options are the following.

=over 4

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are retained. The default
is C<1000>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are retained. The
default is C<10>.


=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir workDir',
                ['lenFilter=i',    'minimum contig length', { default => 1000 }],
                ['covgFilter=f',   'minimum contig mean coverage', { default => 10}],
        );
# Turn off buffering for stdout.
$| = 1;
# Create the loader object.
my $loader = Loader->new();
# Get the statistics object.
my $stats = $loader->stats;
# Save the filter options.
my $lenFilter = $opt->lenfilter;
my $covgFilter = $opt->covgfilter;
# Check the file names.
my ($sampleDir, $workDir) = @ARGV;
if (! $sampleDir) {
    die "Sample directory name missing.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
my ($contigFile, $vectorFile) =
        map { "$sampleDir/$_" } qw(contigs.fasta output.contigs2reads.txt);
if (! -f $contigFile) {
    die "Contig file $contigFile not found.";
} elsif (! -f $vectorFile) {
    die "Vector file $vectorFile not found.";
}
# Check the working directory.
if (! $workDir) {
    $workDir = $sampleDir;
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.";
}
# Now loop through the contig input file, filtering by length.
my %contigs;
print "Processing FASTA file.\n";
my $fh = $loader->OpenFasta(contig => $contigFile);
my $fields = $loader->GetLine(contig => $fh);
while (defined $fields) {
    my ($contig, undef, $seq) = @$fields;
    # Get the sequence length.
    my $len = length $seq;
    # Is this sequence long enough to be meaningful?
    if ($len < $lenFilter) {
        $stats->Add(contigTooShort => 1);
        $contigs{$contig} = 0;
    } else {
        $stats->Add(contigLongEnough => 1);
        $contigs{$contig} = 1;
    }
    # Get the next contig.
    $fields = $loader->GetLine(contig => $fh);
}
# Now compute the coverage information from the vector file.
print "Reading coverage vector file.\n";
my $vh = $loader->OpenFile(coverage => $vectorFile);
# Throw away the first line. It contains headers.
$fields = $loader->GetLine(coverage => $vh);
# Loop through the rest of the file. There will be a mix of nodes and other lines. We only keep the nodes, which
# have contig IDs matching what was in the contig FASTA.
while (! eof $vh) {
    $fields = $loader->GetLine(coverage => $vh);
    my ($contigID, @coverages) = @$fields;
    # Get this contig's status.
    my $status = $contigs{$contigID};
    if (! defined $status) {
        $stats->Add(coverageLineSkipped => 1);
    } else {
        # Compute the mean coverage.
        my $coverage = 0;
        my $count = map { $coverage += $_ } @coverages;
        $coverage /= $count;
        if ($coverage >= $covgFilter) {
            $stats->Add(contigCoverageGood => 1);
            if ($status) {
                $stats->Add(contigGood => 1);
            } else {
                $stats->Add(contigLengthReject => 1);
            }
        } else {
            $stats->Add(contigCoverageBad => 1);
            if ($status) {
                $stats->Add(contigCoverageReject => 1);
            } else {
                $stats->Add(contigBothReject => 1);
            }
        }
    }
}
# Close the vector input file.
close $vh;
# All done.
print "All done.\n" . $stats->Show();
