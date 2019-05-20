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

=head1 Display the Stats Epilog of a Standard Log File

    logtail.pl [ options ]

This script searches for the C<All done.> marker in a log file and displays it along with everything following.  This
usually shows the full statistics dump.

=head2 Parameters

There are no positional parameters.

The command-line options in L<ScriptUtils/ih_options> may be used to modify the standard input.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        );
# Get the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through it, searching for the all-done marker.
my $line = <$ih>;
while (defined $line && $line !~ /^All\s+done/) {
    $line = <$ih>;
}
print "\n";
while (! eof $ih) {
    $line = <$ih>;
    print $line;
}
