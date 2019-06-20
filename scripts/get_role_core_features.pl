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
use RoleParse;

=head1 Get CoreSEED Features for Roles

    get_role_core_features.pl [ options ]

This script will take as input a file of role names and output a file that displays the core feature IDs for each role.
The output file will be a 3-column tab-delimited table containing (0) the role name, (1) the feature ID, and (2) the
name of the genome from which the feature is taken. Each input record could have as many as 1000 output records.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> and L<ScriptUtils::ih_options>.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('fileName', Shrub::script_options(), ScriptUtils::ih_options());
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    # If the input is multi-column, pick just the first column.
    my ($role) = split /\t/, $line;
    # Compute the checksum.
    my $checksum = RoleParse::Checksum($role);
    # Find the role and get all its features.  We use security level 2 to insure we get core-only.
    my @rows = $shrub->GetAll("Role Role2Function Function2Feature Feature2Genome Genome",
                              'Role(checksum) = ? AND Function2Feature(security) = ?',
                              [$checksum, 2], [qw(Function2Feature(to-link) Genome(name))]);
    # Output the results for this role.
    for my $row (@rows) {
        print join("\t", $role, @$row) . "\n";
    }
}