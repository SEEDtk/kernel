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
use gjoseqlib;
use SeedUtils;
use Stats;

=head1 Check DNA Fasta Dumps of Shrub Features

    check_fasta [ options ] file1 file2 ...

This script compares a FASTA file of Shrub database feature DNA to the database protein sequences.
Mismatches will be displayed in the output.

=head2 Parameters

The positional parameters are the names of the FASTA files to be examined. Each file must contain DNA
sequences and the sequence IDs must correspond to Shrub feature IDs.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('file1 file2 ...', Shrub::script_options(),
        ScriptUtils::ih_options(),
        );
my $shrub = Shrub->new_for_script($opt);
my $stats = Stats->new();
for my $file (@ARGV) {
    print "Checking $file.\n-------------------------------\n\n";
    $stats->Add(fileChecked => 1);
    my @triples = gjoseqlib::read_fasta($file);
    for my $triple (@triples) {
        $stats->Add(fidChecked => 1);
        my ($fid, undef, $dna) = @$triple;
        my ($pseq) = $shrub->GetFlat('Feature2Protein Protein', 'Feature2Protein(from-link) = ?', [$fid], 'Protein(sequence)');
        if (! $pseq) {
            $stats->Add(fidNotFound => 1);
        } else {
            my $tseq = SeedUtils::translate($dna, undef, 1);
            $tseq =~ s/\*$//;
            if ($pseq ne $tseq) {
                print "Mismatch for $fid\nP = $pseq\nT = $tseq\n\n";
                $stats->Add(fidMismatch => 1);
            } else {
                $stats->Add(fidMatch => 1);
            }
        }
    }
}
print "All done.\n" . $stats->Show();