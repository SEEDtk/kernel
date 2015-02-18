#!/usr/bin/env perl

#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
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

=head1 Display Roles with EC Numbers

    ShowEcRoles [options] parms

This is a simple script that displays the roles associated with EC numbers. It is used
to verify that the current configuration is working and connects to the database correctly.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options>.

=cut

    use strict;
    use warnings;
    use Shrub;
    use ScriptUtils;

    # Get the command-line options.
    my $opt = ScriptUtils::Opts('', Shrub::script_options());
    # Connect to the database.
    my $shrub = Shrub->new_for_script($opt);
    # Get all the roles with EC numbers and write them out.
    $shrub->PutAll(\*STDOUT, 'Role', 'Role(ec-number) >= ?', ['1'], 'ec-number id description');
