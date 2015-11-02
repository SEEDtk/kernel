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
use Bin;
use Bin::Analyze;
use Shrub;

=head1 Report About Bins

    bins_report.pl [ options ]

Produce a report of the quality bins in a bin list.

=head2 Parameters

The input file is a list of bins,  in JSON format.

The command-line options are those found in L<Shrub::script_options> and L<ScriptUtils/ih_options> plus the following.

=over 4

=item unis

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        Shrub::script_options(),
        ['unis=s',     'universal role file', { default => "$FIG_Config::data/Global/uni_roles.tbl" }],
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
my $binList;
# Read in the universal roles.
my $uniRoles = Bin::Score::ReadUniHash($opt->unifile);
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the analysis object.
my $analyzer = Bin::Analyze->new($shrub);
# This will contain a list of bins for the final statistical report.
my @bins;
# Loop through the bins, collecting them in memory.
while (! eof $ih) {
    my $bin = Bin->new_from_json($ih);
    push @bins, $bin;
}
# Output the report.
$analyzer->Report(\*STDOUT, $uniRoles, \@bins);
