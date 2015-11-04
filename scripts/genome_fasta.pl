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

Display the name of a genome's FASTA file.

=head2 Parameters

The single positional parameter is the genome ID.

The command-line options are those found in L<Shrub/script_options>.
=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('genome_id',
                Shrub::script_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
my $fasta = $shrub->genome_fasta($ARGV[0]);
print "$fasta\n";