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
use Stats;


=head1 Create FASTA for Universal Proteins

    unit_prots.pl [ options ]

This script will output all proteins containing universal functions that are in the L<Shrub> database.  The output FASTA file will
contain the protein ID, the function, and the amino acid sequence.

=head2 Parameters

There are no positional parameters.

The command-line options are those in L<Shrub/script_options> (to select the L<Shrub> database).

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        );
# Create a statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Create the main query.
print STDERR "Reading proteins.\n";
my $q = $shrub->Get("Function Function2Feature Feature2Protein Protein", "Function(universal) = ?", [1], "Function(description) Protein(id) Protein(sequence)");
# This hash will prevent duplicate proteins.
my %pHash;
while (my $record = $q->Fetch()) {
    my ($role, $id, $seq) = $record->Values("Function(description) Protein(id) Protein(sequence)");
    my $count = $stats->Add(protFound => 1);
    if ($pHash{$id}) {
        $stats->Add(protDuplicate => 1);
    } else {
        $pHash{$id} = 1;
        print ">$id $role\n$seq\n";
        $stats->Add(protOut => 1);
    }
    if ($count % 500 == 0) {
        print STDERR "$count proteins read.\n";
    }
}
print STDERR "All done.\n" . $stats->Show();
