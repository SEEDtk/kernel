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

=head1 Display Quality-Related Genome Statistics

    quality_stats.pl [ options ]

This is a pipeline script that adds quality-related information from the database to a tab-delimited file
containing genome IDs.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item col

The index (1-based) of the input column containing the genome IDs. The default is C<0>, indicating the
last column.

=item altFile

The name of a file containing alternative scores. The file must be tab-delimited with genome IDs in the first
column,

=item altColumns

A comma-delimited list of the 1-based indexes of the columns to be copied from the alternative score file.

=back

=head2 Output

The output file is tab-delimited, containing the lines from the input file with the following
additional fields.

=over 4

=item domain

The domain of the organism.

=item dna-size

The number of base pairs in the genome.

=item pegs

The number of protein-encoding genes.

=item well-behaved

C<1> if the genome is well-behaved, else C<0>.

=item core

C<1> if the genome is core, else C<0>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['col|c=i', 'index (1-based) of the input column', { default => 0 }],
        ['altFile=s', 'name of a file containing alternative scores'],
        ['altColumns=s', 'list of alternative-file columns to include']
        );
# Check for an alternative file.
my %altScores;
if ($opt->altfile) {
    my $altFile = $opt->altfile;
    if (! -s $altFile) {
        die "Alternative score file $altFile missing or empty.";
    } elsif (! $opt->altcolumns) {
        die "altColumns required if altFile specified."
    } else {
        my @altColumns = map { $_ - 1 } split /,/, $opt->altcolumns;
        # Read in the alternative-score file.
        open(my $ah, '<', $altFile) || die "Could not open alternative score file: $!";
        while (! eof $ah) {
            my $line = <$ah>;
            chomp $line;
            my @cols = split /\t/, $line;
            $altScores{$cols[0]} = [@cols[@altColumns]];
        }
    }
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the input column index.
my $col = $opt->col;
# Loop through the input file.
while (my @couplets = ScriptUtils::get_couplets($ih, $col, 50)) {
    # Get this batch's genome IDs.
    my @genomes = map { $_->[0] } @couplets;
    # Get the data for these genomes.
    my $filter = 'Genome(id) IN (' . join(', ', map { '?' } @genomes) . ')';
    my %gHash = map { $_->[0] => [@{$_}[1 .. 4]] } $shrub->GetAll('Genome', $filter, \@genomes, 'id domain dna-size well-behaved core');
    # Form them into output.
    for my $couplet (@couplets) {
        my ($genome, $line) = @$couplet;
        my $fields = $gHash{$genome};
        if ($fields) {
            my $alts = $altScores{$genome} // [];
            print join("\t", @$line, @$fields, @$alts) . "\n";
        }
    }
}
