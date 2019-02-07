=head1 Compare Italian Completeness Criterion to EvalG/EvalCon Criterion

    p3x-italian-compare.pl [options]

This script processes the output from L<p3x-process-downloadable-genomes.pl> when applied to the data from the Italian Cell
paper and computes their notion of genome quality (C<L>, C<M>, or C<H>).  It then compares this to our own good-genome flag
in the same output file.

This is an extremely special-purpose script.

Useful statistics will be written to the standard error output.

=head2 Parameters

There are no positional parameters

The standard input can be overridden using the options in L<P3Utils/ih_options> and should be the output file from the
L<p3x-process-downloadable-genomes.pl> script.  The good-genome flag is in a column labeled C<Good?>, and the completeness
and contamination values in columns labeled C<completeness> and C<contamination>.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;

# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::ih_options(),
        ['nohead', 'input file has no headers'],
        );
my $stats = Stats->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.  Note we do not have a key column.
my ($outHeaders) = P3Utils::process_headers($ih, $opt, 1);
# Find the columns of interest.
my (undef, $cols) = P3Utils::find_headers($outHeaders, input => 'completeness', 'contamination', 'Good?');
# Write the output headers.
push @$outHeaders, 'Alt';
P3Utils::print_cols($outHeaders);
# Loop through the input.
while (! eof $ih) {
    # Get the next input line.
    my $line = <$ih>;
    my @fields = P3Utils::get_fields($line);
    my ($complete, $contam, $good) = P3Utils::get_cols(\@fields, $cols);
    $stats->Add(lineIn => 1);
    my $altQuality = 'L';
    if ($contam <= 5.0) {
        if ($complete >= 90) {
            $altQuality = 'H';
        } elsif ($complete >= 50) {
            $altQuality = 'M';
        }
    }
    if ($good) {
        $stats->Add("good-$altQuality" => 1);
        $stats->Add(good => 1);
    } else {
        $stats->Add("bad-$altQuality" => 1);
        $stats->Add(bad => 1);
    }
    $stats->Add("qual-$altQuality" => 1);
    P3Utils::print_cols([@fields, $altQuality]);
}
print STDERR "All done.\n" . $stats->Show();
