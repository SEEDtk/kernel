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
use Shrub::Functions;
use RoleParse;

=head1 Validate Roles

    role_check.pl [ options ]

This script reads a table of functional assignments and validates the roles found.
It will output the rows containing assignments whose roles are malformed or not found in the database.

=head2 Parameters

There are no positional parameters.

The input file is tab-delimited, with functional assignments in the last column.

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item col

The column number (1-based) containing the functional assignments. The default is C<0>, which returns the
last column.

=item hypo

Include hypothetical proteins in the output.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ["col|c=i", 'column containing functional assignment text'],
        ["hypo|h", 'if specified, hypothetical proteins will be included'],
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# We need to loop through the input, filtering out rows with good roles in them.
while (! eof $ih) {
    my @couplets = ScriptUtils::get_couplets($ih, $opt->col, 10);
    for my $couplet (@couplets) {
        my ($func, $row) = @$couplet;
        # Process this function to get the list of roles.
        my ($statement, undef, $roles) = Shrub::Functions::Parse($func);
        # This will be set to TRUE if we can't find a role.
        my $badLine;
        if ($statement =~ /hypothetical/) {
            # Hypothetical roles are ok. We don't bother to process them, though,
            # unless the user wants to see them.
            $badLine = $opt->hypo;
        } elsif (! @$roles) {
            # Here the function is malformed.
            $badLine = 1;
        } else {
            # Here we need to check the roles.
            for my $role (@$roles) {
                my $checksum = RoleParse::Checksum($role);
                # Compute the ID for this checksum.
                my ($id) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checksum], 'id');
                if (! $id) {
                    $badLine = 1;
                }
            }
        }
        # If we are malformed or found a bad role, output the line.
        if ($badLine) {
            print join("\t", @$row), "\n";
        }
    }
}
