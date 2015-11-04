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

=head1 Display Fasta File Name

    genome_fasta.pl [ options ] genome_id

Display the name or content of a genome's FASTA file.

=head2 Parameters

The single positional parameter is the genome ID.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item pipe

If specified, the FASTA data will be written to the standard output instead of the file name.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('genome_id',
                Shrub::script_options(),
                ['pipe|p', 'write FASTA data to standard output']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the genome ID.
my ($genome) = @ARGV;
if (! $genome) {
    die "A genome ID is required.";
}
# Get the FASTA file name.
my $fasta = $shrub->genome_fasta($genome);
if (! $fasta) {
    die "Genome $genome not found.";
}
if ($opt->pipe) {
    # Here we are writing the FASTA data to the output stream.
    open(my $ih, "<$fasta") || die "Could not open fasta file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        print $line;
    }
} else {
    # Here we just want the name.
    print "$fasta\n";
}
