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


package ScriptUtils;

    use strict;
    use warnings;

=head1 Script Utilities

This package contains utilities that are useful in coding SEEDtk command-line scripts.

=head2 Public Methods

=head3 IH

    my $ih = ScriptUtils::IH($fileName);

Get the input file handle. If the parameter is undefined or empty, the
standard input will be used. Otherwise the file will be opened and an
error thrown if the open fails. When debugging in Eclipse, the
standard input is not available, so this method provides a cheap way for
the input to ber overridden from the command line. This method provides a
compact way of insuring this is possible. For example, if the script has
two positional parameters, and the last is an optional filename, you
would code

    my $ih = ScriptUtils::IH($ARGV[1]);

If the C<-i> option contains the input file name, you would code

    my $ih = ScriptUtils::IH($opt->i);

=over 4

=item fileName

Name of the file to open for input. If the name is empty or omitted, the
standard input will be returned.

=item RETURN

Returns an open file handle for the script input.

=back

=cut

sub IH {
    # Get the parameters.
    my ($fileName) = @_;
    # Declare the return variable.
    my $retVal;
    if (! $fileName) {
        # Here we have the standard input.
        $retVal = \*STDIN;
    } else {
        # Here we have a real file name.
        open($retVal, "<$fileName") ||
            die "Could not open input file $fileName: $!";
    }
    # Return the open handle.
    return $retVal;
}


1;
