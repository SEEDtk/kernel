=head1 Analyze Drug Groups by Various Criteria

    get_drug_groups.pl [options] criterion

This program will isolate the input file into drug groups. Each drug group represents a pair of drugs applied to a cell line. We then
apply the specified criterion to the drug group's elements to assign a score. The N groups with the highest score are output.

=head2 Parameters

The positional parameter is the name of the criterion to apply. Currently, the following criteria are supported.

=over 4

=item deltaDiff

The score is the largest difference between growth rates for adjacent drug concentrations for one of the drugs (that is, the other drug's
concentration is constant). The difference is computed as the percent difference relative to the larger growth rate.

=item synergy

The score is the largest difference between the control growth rate and the actual growth rate.

=item synergy

This is the same as synergy, but rows with negative actual growth rates are discarded.

=item antagony

The score is the smallest different between the control growth rate and the actual growth rate.

=item antagony2

This is the same as antagony, but rows with negative actual growth rates are discarded.

=back

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The following additional command-line options are supported.

=over 4

=item N

Number of results to return.

=back

=cut

use strict;
use P3Utils;

# Get the command-line options.
my $opt = P3Utils::script_opts('criterion', P3Utils::ih_options(),
        ['N=i', 'number to return', { default => 10 }],
        );
my ($criterion) = @ARGV;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Process the header.
my $line = <$ih>;
$line =~ s/[\r\n]+$//;
my @headers = split /,/, $line;
# This array tracks the best N groups. It is sorted by score; each tuple is a score followed by a string for the group.
my @saved;
P3Utils::print_cols(['score', @headers]);
# This is a counter for tracing.
my $count = 0;
# This array contains the lines for the current group.
my $lines;
# This contains the name of the current group.
my $groupName = '';
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    my @cols = split /,/, $line;
    my $group = "$cols[8],$cols[14],$cols[28]";
    if ($group ne $groupName) {
        ProcessGroup($groupName, $lines, \@saved, $criterion) if $groupName;
        $groupName = $group;
        $lines = [];
    }
    push @$lines, \@cols;
}
# We don't process the residual-- the file comes with a trailer line.
# Print the saved groups.
for my $groupData (@saved) {
    my ($score, $group) = @$groupData;
    for my $line (@$group) {
        print join("\t", $score, @$line) . "\n";
    }
}

# Score a group and merge it into the saved list.
sub ProcessGroup {
    my ($groupName, $lines, $saved, $criterion) = @_;
    if (++$count % 1000 == 0) {
        print STDERR "Processing group $count: $groupName.\n";
    }
    # Compute the score.
    my $score;
    if ($criterion eq 'deltaDiff') {
        $score = DeltaDiff($lines);
    } elsif ($criterion eq 'synergy') {
        $score = Synergy($lines);
    } elsif ($criterion eq 'synergy2') {
        $score = Synergy2($lines);
    } elsif ($criterion eq 'antagony') {
        $score = Antagony($lines);
    } elsif ($criterion eq 'antagony2') {
        $score = Antagony2($lines);
    } else {
        die "Unknown criterion $criterion.";
    }
    # Only proceed if we have a score.
    if (defined $score) {
        # Sort the new group into the saved set.
        @$saved = sort { $b->[0] <=> $a->[0] } (@$saved, [$score, $lines]);
        # Insure we are not too big.
        if (scalar @$saved > $opt->n) {
            pop @$saved;
        }
    }
}

# Useful columns when scoring.
use constant { CONCIDX1 => 10, CONCIDX2 => 16, GROWTH => 19, EXPECTED => 24 };

# Score the group based on percent difference between adjacent concentrations.
# The group is a list of lists. Each sub-list is a line.
sub DeltaDiff {
    my ($lines) = @_;
    # Create a matrix of growth rates.
    my @growths;
    my $maxIdx = 0;
    for my $line (@$lines) {
        # Only proceed if we have a good growth rate.
        if ($line->[EXPECTED()] ne '') {
            my ($idx1, $idx2) = ($line->[CONCIDX1()] - 1, $line->[CONCIDX2()] - 1);
            $maxIdx = $idx1 if $idx1 > $maxIdx;
            if ($idx2 < 0) {
                # Here we have a single drug, not a combo.
                $idx2 = 0;
            }
            $maxIdx = $idx2 if $idx2 > $maxIdx;
            $growths[$idx1][$idx2] = $line->[EXPECTED()];
        }
    }
    # Now we have a two-dimensional array of growth rates, ordered by concentration in each dimension. Find the
    # largest percent difference. We do this once in each direction.
    my $retVal = 0;
    for (my $idx1 = 0; $idx1 <= $maxIdx; $idx1++) {
        for (my $idx2 = 1; $idx2 <= $maxIdx; $idx2++) {
            my $score = pct($growths[$idx1][$idx2 - 1], $growths[$idx1][$idx2]);
            $retVal = merge($score, $retVal);
        }
    }
    for (my $idx2 = 0; $idx2 <= $maxIdx; $idx2++) {
        for (my $idx1 = 1; $idx1 <= $maxIdx; $idx1++) {
            my $score = pct($growths[$idx1-1][$idx2], $growths[$idx1][$idx2]);
            $retVal = merge($score, $retVal);
        }
    }
    # Return the computed score.
    return $retVal;
}

# Synergy score. Expected - Actual.
sub Synergy {
    my ($lines) = @_;
    my $retVal;
    for my $line (@$lines) {
        my ($expected, $actual) = ($line->[EXPECTED()], $line->[GROWTH()]);
        if (defined $expected && defined $actual) {
            my $score = $expected - $actual;
            $retVal = merge($score, $retVal);
        }
    }
    return $retVal;
}

# Modified Synergy score. Expected - Actual, Actual >= 0
sub Synergy2 {
    my ($lines) = @_;
    my $retVal;
    for my $line (@$lines) {
        my ($expected, $actual) = ($line->[EXPECTED()], $line->[GROWTH()]);
        if (defined $expected && defined $actual && $actual >= 0) {
            my $score = $expected - $actual;
            $retVal = merge($score, $retVal);
        }
    }
    return $retVal;
}

# Antagonism score. Actual - Expected.
sub Antagony {
    my ($lines) = @_;
    my $retVal;
    for my $line (@$lines) {
        my ($expected, $actual) = ($line->[EXPECTED()], $line->[GROWTH()]);
        if (defined $expected && defined $actual) {
            my $score = $actual - $expected;
            $retVal = merge($score, $retVal);
        }
    }
    return $retVal;
}

# Modified antagonism score. Actual - Expected, Actual >= 0
sub Antagony2 {
    my ($lines) = @_;
    my $retVal;
    for my $line (@$lines) {
        my ($expected, $actual) = ($line->[EXPECTED()], $line->[GROWTH()]);
        if (defined $expected && defined $actual && $actual >= 0) {
            my $score = $actual - $expected;
            $retVal = merge($score, $retVal);
        }
    }
    return $retVal;
}

# Balanced percent difference. Return undef if data is missing.
sub pct {
    my ($a, $b) = @_;
    my ($retVal, $base);
    if (defined $a && defined $b) {
        if (abs($a) < abs($b)) {
            $base = $b;
        } else {
            $base = $a;
        }
        if ($base > 0) {
            $retVal = abs($a - $b) * 100 / $base;
        }
    }
    return $retVal;
}

# Merge a new score into the best value.
sub merge {
    my ($score, $best) = @_;
    my $retVal = $best;
    if (defined $score && (! defined $best || $score > $best)) {
        $retVal = $score;
    }
    return $retVal;
}
