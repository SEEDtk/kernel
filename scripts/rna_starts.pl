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
use BasicLocation;
use FastA;

=head1 Produce a Training File for RNA Starts

    rna_starts.pl [ options ] fastaFile tblFile

This script creates a training set for RNA start detection.  The positional parameters are the names of a FASTA file
containing proteins and a flat file containing BLAST matches.  Both files should be the output of a <genome.distance match>
command.

The FASTA file contains the protein sequences found in the RNA along with the upstream DNA for each.  The label for each
sequence is of the form I<sequenceID>C<.p>I<pos>C<.F>I<frame>.  In the flat file, the first column is the label and the
fourth column is the sequence length in codons.  We remember the label and the full extent of each protein in the sequence.
For the FASTA file, a sequence with a label found in the flat file is considered a hit.  A sequence that starts inside one
of the hit sequences is discarded.  All other sequences are misses.  We then use the upstream DNA to form the training set.

=head2 Parameters

The positional parameters are the name of the FASTA file and the name of the flat file.  The training set is written to
standard output and progress messages are written to the standard error log.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('fastaFile tblFile',
        );
# Create a statistics object.
my $stats = Stats->new();
# Get the files.
my ($fastaFile, $tblFile) = @ARGV;
if (! $fastaFile) {
    die "No FASTA file specified.";
} elsif (! $tblFile) {
    die "No flat file specified.";
}
# This will store the labels for the hits.
my %hits;
# For each contig, this will list the sequence locations.
my %locs;
# Loop through the flat file.
open(my $ih, '<', $tblFile) || die "Could not open $tblFile: $!";
# Skip the header line.
my $line = <$ih>;
while (! eof $ih) {
    $line = <$ih>;
    $stats->Add(lineIn => 1);
    if ($line =~ /^(\S+)\.p(\d+)\.F(\d+)\t\S+\t[^\t]+\t(\d+)\t/) {
        my ($contig, $start, $frm, $len) = ($1, $2, $3, $4);
        $stats->Add(hitLine => 1);
        $hits{"$contig.p$start.F$frm"} = 1;
        push @{$locs{$contig}}, BasicLocation->new($contig, $start, '+', ($len+1)*3);
    }
}
# Close the flat file.
close $ih; undef $ih;
# Write the output header.
print join("\t", 'label', 'hit', map { "p$_" } 1 .. 50) . "\n";
# The output lines will be stored in here.
my %lines;
# Open the FASTA file.
my $fh = FastA->new($fastaFile);
while ($fh->next) {
    $stats->Add(seqIn => 1);
    # Analyze this sequence.
    my $label = $fh->id();
    if ($hits{$label}) {
        $stats->Add(hitFound => 1);
        WriteHit("coding", $fh);
    } else {
        my ($contig, $pos) = ($label =~ /^(\S+)\.p(\d+)\./);
        my @locs = @{$locs{$contig}};
        my $ok = 1;
        while (scalar @locs > 0 && $ok) {
            my $loc = pop @locs;
            if ($loc->Left <= $pos && $loc->Right >= $pos) {
                $ok = 0;
            }
        }
        if (! $ok) {
            $stats->Add(discarded => 1);
        } else {
            $stats->Add(missFound => 1);
            WriteHit("other", $fh);
        }
    }
}
# Spool off the lines.
my $hitLines = $lines{"coding"};
my $missLines = $lines{"other"};
while (scalar @$hitLines) {
    print pop @$hitLines;
    if (scalar @$missLines) {
        print pop @$missLines;
    }
    if (scalar @$missLines) {
        print pop @$missLines;
    }
}
while (scalar @$missLines) {
    print pop @$missLines;
}
print STDERR "All done.\n" . $stats->Show();

# Output a hit.
sub WriteHit {
    my ($type, $fh) = @_;
    push @{$lines{$type}}, join("\t", $fh->id, $type, split //, $fh->comment) . "\n";
    $stats->Add(lineOut => 1);
}
