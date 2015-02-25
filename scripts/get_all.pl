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
use Shrub;
use ScriptUtils;

=head1 Perform a General Database Query

    get_all [ options ] -p path -c ' constraint' -v parm1 -v parm2 ... field1 field2 ...

This script performs a general database query. It can either perform a standalone query
or process data from a tab-delimited input file. If the input is from a tab-delimited
input file, the output file will consist of the result columns tacked on to the
end of each input row.

The L<ERDB/GetAll> command is used to perform the query processing, and the
documentation there should be consulted for more details.

Essentially, to navigate the database you must specify a I<path> through the entities and
relationships, a I<constraint> that limits the output, one or more I<parameters> to the
constraint, and a list of I<fields> to be returned.

=head2 Parameters

The positional parameters are the names of the fields to be returned in order.
The field names should be specified in L<ERDB/Standard Field Name Format>.

The command-line options are those found in L<Shrub/script_options>
and L<ScriptUtils/ih_options> plus the following.

=over 4

=item path

A list of entity and relationship names, space-delimited, describing a
path through the database. If there is more than one name in the path,
you will have to enclose it in quotes.

=item constraint

A query constraint, in the form of an SQL WHERE clause (without the word
C<WHERE>). The database field names must be specified in L<ERDB/Standard Field Name Format>,
and you can specify parameterizations using question marks. The question marks
will be replaced by the values of the I<value> option, in order.

=item value

A list of the values to be substituted for the parameter marks. This is an option
array, so if you need to specify more than one value, you simply specify the
option more than once. If you are using an input file, you can specify a value
taken from a column of the current input record using a dollar sign (C<$>)
followed by a column number. The column numbers are 1-based so a value parameter of
C<$1> would indicate the first column of the input file. The special value C<$n>
indicates the last column of the input file.

=item count

Maximum number of records to return. The default is C<0>, meaning all records
will be returned. If an input file is being used, this constraint applies to
each row of input. So, a count of 5 with an input file of 10 records would mean
a maximum of 50 results.

=back

=head2 Examples

The following command reads a list of subsystem IDs from the standard input
and lists all the genomes in each subsystem along with each genome's variant
code. The subsystem IDs are taken from the first column of the input.

    get_all -p "Subsystem Subsystem2Row Row2Genome" -c "Subsystem(id) = ?" \
            -v \$1 "Row2Genome(to-link)" "Row2Genome(variant)"

The following command does the same thing, but instead of reading from the
standard input it just displays the genomes for Histidine Degradation.
Note that we include the Subsystem ID in the output so that we have the
same three columns as in the previous example.

    get_all -p "Subsystem Subsystem2Row Row2Genome" -c "Subsystem(id) = ?" \
            -v "HistDegr" "Subsystem(id)" "Row2Genome(to-link)" \
            "Row2Genome(variant)"


This reads a list of protein MD5s from the standard input and outputs
the IDs of the features that produce the protein. The protein MD5s
are in the last column of the input file.

    get_all -p "Protein2Feature" -c "Protein2Feature(from-link) = ?" \
            -v \$n "Protein2Feature(to-link)"

=cut

# Get the command parameters.
my $opt = ScriptUtils::Opts('field1 field2 ...', Shrub::script_options(),
        ScriptUtils::ih_options(),
        ["path|p=s", "path through the database", { required => 1} ],
        ["constraint|c=s", "query constraint (if any)"],
        ["value|v=s@", "parameter values for the query constraint (multiple)"],
        ["count|n=i", "maximum number of records to output per input line (default is 0, meaning return everything)", { default => 0 }]
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the list of parameter values. If none were supplied, we use an empty list.
my $valueList = $opt->value // [];
# Check for input markers. This list will contain undef if the parameter value is a constant,
# and a column index if the parameter value is taken from the input line.
my @inputList;
# This will be set to the input file handle if we find at least one input marker.
my $ih;
# Loop through the parameter values.
my $n = scalar @$valueList;
for (my $i = 0; $i < $n; $i++) {
    if ($valueList->[$i] =~ /^\$(\d+|n)/) {
        # We have an input column marker. Save the column number.
        my $col = (($1 eq 'n') ? -1 : ($1 - 1));
        push @inputList, ($1 - 1);
        # If this is our first marker, open the input file.
        if (! defined $ih) {
            $ih = ScriptUtils::IH($opt->input);
        }
    } else {
        # Here we have an ordinary constant value.
        push @inputList, undef;
    }
}
# Now get the path, constraint, and count. The count already defaults to 0 so we don't
# need to worry about it being undefined.
my $path = $opt->path;
my $constraint = $opt->constraint // '';
my $count = $opt->count;
# Get the list of output field names.
my @fields = @ARGV;
if (! scalar(@fields)) {
    die "No output fields specified.\n";
}
# Get the first line of input. If we have no input, this is an empty string.
my $line = (defined $ih ? <$ih> : "\n");
# Loop until we run out of input.
while (defined $line) {
    # Parse the input line.
    chomp $line;
    my @cols = split /\t/, $line;
    # Form the list of parameter values.
    my @parms;
    for (my $i = 0; $i < $n; $i++) {
        if (defined $inputList[$i]) {
            my $colIdx = $inputList[$i];
            $colIdx = $#cols if $colIdx < 0;
            push @parms, $cols[$colIdx];
        } else {
            push @parms, $valueList->[$i];
        }
    }
    # Query the database.
    my @rows = $shrub->GetAll($path, $constraint, \@parms, \@fields, $count);
    # Output the results.
    for my $row (@rows) {
        print join("\t", @cols, @$row) . "\n";
    }
    # Get the next input line (or undef if there is no input file).
    $line = (defined $ih ? <$ih> : undef);
}
# All done. Close the input.
if (defined $ih) {
    close $ih;
}
