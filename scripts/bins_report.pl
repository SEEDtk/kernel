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

=head1 Compute Statistics About Bins

    bins_report.pl [ options ]

Produce statistics about a set of bins.

=head2 Parameters

The input file is a list of bins, either in L<bin exchange format|Bin/Bin Exchange Format> or in JSON format.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item exchange

If specified, it is presumed the bins are in exchange format. Otherwise, they are read in JSON format.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        ['exchange|x', 'exchange format input']
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
my $binList;
# Read the bins in the specified format.
if ($opt->exchange) {
    $binList = Bin::ReadContigs($ih);
} else {
    $binList = Bin::ReadBins($ih);
}
# Get the bins report.
my $stats = Bin::Analyze::Report($binList);
print $stats->Show();
