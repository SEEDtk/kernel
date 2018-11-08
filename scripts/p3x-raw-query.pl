=head1 Query Test Script

    p3x-raw-query.pl [options] table id-name attr1 attr2 ... attrN

This script performs a raw P3 query that can be used to analyze the contents of the SOLR database.  It is for development
purposes only.

=head2 Parameters

The positional parameters are the name of the table to query followed by the table's ID field and then a list of the other
fields to be output (this list can be empty).

The following command-line options are supported.

=over 4

=item limit

The number of records to return.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;


# Get the command-line options.
my $opt = P3Utils::script_opts('table id-name attr1 attr2 ... attrN',
        ['limit=i', 'number of records to return', { default => 100 }]
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Parse the parameters.
my ($table, @fields) = @ARGV;
# Write the header.
P3Utils::print_cols(\@fields);
# Form the query.
my @results = $p3->query($table, ['ne', $fields[0], 'X'], ['select' => @fields], [limit => $opt->limit]);
# Process the results.
for my $result (@results) {
    my @values = map { $result->{$_} } @fields;
    P3Utils::print_cols(\@values);
}
