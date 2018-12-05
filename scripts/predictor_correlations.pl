=head1 Compute Correlations in a Predictor Directory

    predictor_correlations.pl [options] predDir

This script reads a predictor directory (C<row.h>, C<col.h>, C<X>) and outputs simple percentage correlations.
A column object's correlation with a second column object is the number of times the second object occurs in
a row where the first occurs divided by the number of rows in which the first occurs.  This is a very simple
computation to make from a predictor directory format.

=head2 Parameters

The single positional parameter is the name of the input predictor directory.

Additional command-line options are as follows.

=over 4

=item verbose

Progress messages will be written to STDERR.

=item min

The minimum correlation percentage required to output a pair.  The default is C<80>.

=item minCount

The minimum number of rows in which a column object must occur for its correlation to count. The default is C<5>.

=back

=cut

use strict;
use P3Utils;
use Math::Round;

# Get the command-line options.
my $opt = P3Utils::script_opts('predDir',
        ['verbose|debug|v', 'show progress on STDERR'],
        ['min|m=f', 'minimum percentage for correlation', { default => 80 }],
        ['minCount|k=i', 'minimum number of occurrences for correlation to count', { default => 5 }]
        );
# Get the input directory.
my ($predDir) = @ARGV;
if (! $predDir) {
    die "No input directory specified.";
} elsif (! -d $predDir) {
    die "$predDir not found or invalid.";
} elsif (! -s "$predDir/X" || ! -s "$predDir/col.h") {
    die "$predDir does not appear to be a predictor directory.";
}
# Get the parameters.  Note we insure the minimum count is at least 1.
my $debug = $opt->verbose;
my $min = $opt->min;
my $minCount = $opt->mincount;
$minCount = 1 if $minCount < 1;
# Read in the column names.  This is tricky because we support multi-field names.
# The following variables will contain the output column headers to use for each object name.
my ($col1, $col2);
# The following contains the name of each column, in order.
my @cols;
# Now we read the file.
print STDERR "Reading column names.\n" if $debug;
my $maxFlds = 0;
open(my $ih, "<$predDir/col.h") || die "Could not open col.h: $!";
while (! eof $ih) {
    my $line = <$ih>;
    my ($idx, @name) = P3Utils::get_fields($line);
    if (scalar(@name) > $maxFlds) {
        $maxFlds = scalar @name;
    }
    $cols[$idx] = join("\t", @name);
}
print STDERR scalar(@cols) . " columns found.\n" if $debug;
close $ih; undef $ih;
# Now create the column labels.
my @suffix;
if ($maxFlds > 1) {
    push @suffix, '' x ($maxFlds - 1);
}
$col1 = join("\t", "Column 1", @suffix);
$col2 = join("\t", "Column 2", @suffix);
# Output the headers.
print join("\t", $col1, $col2, "Percent", "Count") . "\n";
# This hash will contain the correlation counts, keyed by col1 index then col2 index.
my %counts;
# This array tracks how many rows each column occurs in.
my @counts;
print STDERR "Reading matrix.\n" if $debug;
my $rowIdx = 0;
open($ih, "<$predDir/X") || die "Could not open X: $!";
while (! eof $ih) {
    my $line = <$ih>;
    my @cells = P3Utils::get_fields($line);
    # This will record all the columns found in the row.
    my %cols;
    for (my $colIdx = 0; $colIdx < @cells; $colIdx++) {
        my $val = $cells[$colIdx];
        if ($val > 0) {
            $cols{$colIdx} = 1;
        }
    }
    # Now we count all the occurrences.
    for my $colIdx (keys %cols) {
        $counts[$colIdx]++;
        for my $colIdx2 (keys %cols) {
            if ($colIdx2 != $colIdx) {
                $counts{$colIdx}{$colIdx2}++;
            }
        }
    }
    # Update the row counter.
    $rowIdx++;
    print STDERR "$rowIdx rows processed.\n" if $debug && ($rowIdx % 200 == 0);
}
# Now we compute the percentages in here.  The key will be the concatenation of column names.
my %percents;
print STDERR "Calculating percentages.\n" if $debug;
my $computed = 0;
for my $colIdx (keys %counts) {
    $col1 = $cols[$colIdx];
    my $denom = $counts[$colIdx] // 0;
    if ($denom >= $minCount) {
        my $subCounts = $counts{$colIdx};
        for my $colIdx2 (keys %$subCounts) {
            my $percent = $subCounts->{$colIdx2} * 100 / $denom;
            if ($percent >= $min) {
                my $name = "$col1\t$cols[$colIdx2]";
                $percents{$name} = [Math::Round::nearest(0.01, $percent), $subCounts->{$colIdx2}];
            }
            $computed++;
            print STDERR "$computed percentages computed.\n" if $debug && ($computed % 1000 == 0);
        }
    }
}
print STDERR "Sorting results.\n" if $debug;
my @keys = sort { $percents{$b}[0] <=> $percents{$a}[0] || $percents{$b}[1] <=> $percents{$a}[1] } keys %percents;
print STDERR "Printing results.\n" if $debug;
for my $key (@keys) {
    print join("\t", $key, @{$percents{$key}}) . "\n";
}
