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
use KmerDb;

=head1 Output Kmers That Connect Clusters

    cluster_detail_kmers.pl [ options ] kmer_db

This script identifies the kmers in a kmer database that bind together sequence groups. The standard input contains
group IDs, and for each one, the ID of each connected group and the kmer that connects them is output.

=head2 Parameters

The single positional parameter is the file name of the kmer database.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item col

Index (1-based) of the column containing the input group ID. A value of C<0> indicates the last column.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('kmer_db', ScriptUtils::ih_options(),
        ['col|c=i', 'index (1-based) of input column containing group IDs', { default => 0 }],
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
my $col = $opt->col;
# Create a hash of the incoming group IDs. This hash maps each group ID to its input line.
my %inputs;
while (! eof $ih) {
    my ($couplet) = ScriptUtils::get_couplets($ih, $col, 1);
    my ($id, $line) = @$couplet;
    $inputs{$id} = $line;
}
# Create the kmer database.
my $kmerdb;
my ($kmerFile) = @ARGV;
if (! $kmerFile) {
    die 'A kmer database name is required.';
} elsif (! -f $kmerFile) {
    die "Invalid kmer file name: $kmerFile.";
} else {
    $kmerdb = KmerDb->new(json => $kmerFile);
    # Now we loop through the kmers, finding connections.
    my ($kmer, $groupList);
    while (($kmer, $groupList) = each %{$kmerdb->{kmerHash}}) {
        # Get this kmer's groups.
        my @sorted = sort @$groupList;
        # Only proceed if it's at least a double.
        if (scalar(@sorted) > 1) {
            # Look for interesting groups.
            for my $group (@sorted) {
                my $line = $inputs{$group};
                if ($line) {
                    # Here the group is interesting. Output the connections.
                    for my $group2 (@sorted) {
                        # Insure this connection has not already been output.
                        if (! $inputs{$group2} || $group2 gt $group) {
                            print join("\t", @$line, $group2, $kmer) . "\n";
                        }
                    }
                }
            }
        }
    }
}
