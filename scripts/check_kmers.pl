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

=head1 Check Signature Kmers

    check_kmers.pl [ options ] repFile

This script examines signature kmers to determine how many correspond to each representative genome.
It takes as input a signature kmer file and the output from L<rep_server.pl> describing the representative
genomes themselves. The output file will be a copy of the signature kmer file with the signature counts
appended. This file can then be read into L<rep_server.pl>.

=head2 Parameters

The positional parameter is the name of the representative genomes file. The standard input should be the signature
Kmers file. Both files are tab-delimited. The signature kmers file has one record per kmer, containing (0) the kmer
itself and (1) the associated representative genome's ID. The representative genomes file has one record per genome,
containing (0) a sequence number, (1) the genome ID, and (2) the genome name.

The command-line options are those found in L<ScriptUtils/ih_options>, specifying the standard input.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('repFile', ScriptUtils::ih_options(),
        );
my $stats = Stats->new();
# Initialize a hash for the representative genomes.
print STDERR "Processing representatives.\n";
my %repCounts;
my %repNames;
my ($repFile) = @ARGV;
if (! $repFile) {
    die "No representative genomes file specified.";
} elsif (! -s $repFile) {
    die "Representative genomes file $repFile missing or empty.";
} elsif (! open(my $rh, "<$repFile")) {
    die "Could not open $repFile: $!";
} else {
    while (! eof $rh) {
        my $line = <$rh>;
        $stats->Add(repLineIn => 1);
        if ($line =~ /(\d+\.\d+)\t(.+)/) {
            my ($id, $name) = ($1, $2);
            $stats->Add(genomesIn => 1);
            $repNames{$id} = $name;
            $repCounts{$id} = 0;
        }
    }
}
# Now we have a count of zero for each representative, plus a hash to give us the names for
# output purposes. Open the input file. Note we echo to the output.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the input, counting kmers.
while (! eof $ih) {
    my $line = <$ih>;
    print $line;
    $stats->Add(kmerLineIn => 1);
    if ($line =~ /^([acgt]+)\t(\d+\.\d+)/) {
        my ($kmer, $id) = ($1, $2);
        $stats->Add(kmersIn => 1);
        if (! exists $repNames{$id}) {
            $stats->Add(invalidGenomeID => 1);
        } else {
            $repCounts{$id}++;
            my $count = $stats->Add(kmerCounted => 1);
            if ($count % 100000 == 0) {
                print STDERR "$count kmers counted.\n";
            }
        }
    }
}
print STDERR "Sorting kmers.\n";
# Put a trailer in the output.
print "//\n";
# Now we have the number of kmers for each genome. Sort them and write them out.
my @genomes = sort { $repCounts{$b} <=> $repCounts{$a} } keys %repNames;
for my $id (@genomes) {
    if (! $repCounts{$id}) {
        $stats->Add(kmerFreeGenome => 1);
    }
    print "$id\t$repCounts{$id}\t$repNames{$id}\n";
}
print STDERR "All done.\n" . $stats->Show();
