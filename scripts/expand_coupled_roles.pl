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

=head1 Convert Roles to Descriptions in a Coupling File

    expand_coupled_roles.pl [ options ] 

This script reads the output of a coupled-roles script and converts the role IDs to names.

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(), ScriptUtils::ih_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Echo the header.
my $line = <$ih>;
print $line;
# We'll cache role IDs in here.
my %roles;
# Convert the other lines.
while (! eof $ih) {
    my $line = <$ih>;
    # Note the count includes the new-line, but we don't care.
    my ($role1, $role2, $count) = split /\t/, $line;
    $role1 = get_role($role1, \%roles);
    $role2 = get_role($role2, \%roles);
    print join("\t", $role1, $role2, $count);
}

sub get_role {
    my ($role, $rolesH) = @_;
    my $retVal = $rolesH->{$role};
    if (! $retVal) {
        ($retVal) = $shrub->GetFlat('Role', 'Role(id) = ?', [$role], 'description');
        $rolesH->{$role} = $retVal;
    }
    return $retVal;
}
