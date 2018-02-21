=head1 Filter Samples to Remove Poor Quality

    p3-filter-samples.pl [options] <inFile >outFile 2>logFile

This script will take as input a list of SRA sample IDs and filter out the ones unlikely to produce good bins.
The basic strategy is to count the spots and throw away the sample if it has fewer than 1 million.


=head2 Parameters

There are no positional parameters.


The options in L<P3Utils/col_options> can be used to specify the input column and
L<P3Utils/ih_options> can be used to modify the standard input.

The following additional options are supported.

=over 4

=item min

The minimum number of spots required. The default is C<1000000>.

=back

=cut

use strict;
use P3Utils;
use Stats;
use File::Copy::Recursive;
use SRAlib;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::col_options(), P3Utils::ih_options(),
        ['min=i', 'minimum number of spots', { default => 1000000 }],
        );
# Create a statistics object.
my $stats = Stats->new();
# Get the options.
my $min = $opt->min;
# Open the standard input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Write the output headers.
if (! $opt->nohead) {
    P3Utils::print_cols($outHeaders);
}
# Create the SRAlib object.
my $sra = SRAlib->new(logH => \*STDERR, stats => $stats);
# Loop through the samples.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    for my $couplet (@$couplets) {
        my ($id, $line) = @$couplet;
        my ($spots) = $sra->get_stats($id);
        if (! $spots) {
            print STDERR "$id not found.\n";
            $stats->Add(notFound => 1);
        } elsif ($spots < $min) {
            print STDERR "$id too small-- $spots spots.\n";
            $stats->Add(tooSmall => 1);
        } else {
            P3Utils::print_cols($line);
            $stats->Add(accepted => 1);
        }
    }
}
print STDERR "All done:\n" . $stats->Show();
