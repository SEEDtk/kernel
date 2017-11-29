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

The score is the largest arithmetic difference between the control growth rate and the actual growth rate.

=item synergy

This is the same as synergy, but rows with negative actual growth rates are discarded.

=item antagony

The score is the smallest arithmetic difference between the control growth rate and the actual growth rate.

=item antagony2

This is the same as antagony, but rows with negative actual growth rates are discarded.

=back

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The following additional command-line options are supported.

=over 4

=item N

Number of results to return. The default is 10.

=item min

If specified, ALL results with a score at or above the specified value will be returned. In this case, the B<N> option will be ignored.

=item summary

If specified, the name of a file into which the number of results output for each pair of drugs will be written.

=item doubles

If specified, only groups with two drugs will be included in the output.

=item grouped

If specified, the output will be sorted by drug group instead of score.

=item names

If specified, the name of a file mapping drug numbers to names. The file should be tab-delimited, with numbers in the first column and
names in the second, and no headers.

=item significant

The number of results considered significant in the summary report. Score counts lower than this value will not be shown.

=item triples

If specified, the name of a file to contain a report on significant drug triplets.

=item tripsig

The number of scores considered significant for triplets. A set of three drugs is a triplet if each pair in the triplet has this number of
scores or greater.

=item panels

If specified, the name of a file to contain a report on the types of cell lines that are considered significant for each drug.

=back

=cut

use strict;
use P3Utils;

# Get the command-line options.
my $opt = P3Utils::script_opts('criterion', P3Utils::ih_options(),
        ['N=i', 'number to return', { default => 10 }],
        ['min=f', 'return all groups with a score ge this value'],
        ['doubles', 'only include groups with two drugs'],
        ['grouped', 'sort output by drug rather than score'],
        ['summary=s', 'summary file name'],
        ['names=s', 'file mapping drug numbers to names'],
        ['significant=i', 'the number of results considered significant in the summary report', { default => 1 }],
        ['triples=s', 'triples report file name'],
        ['tripsig=i', 'the number of results considered significant in the triples report', { default => 7 }],
        ['panels=s', 'cell lines type report file name'],
        );
my ($criterion) = @ARGV;
# Check for drug names.
my %nameH;
if ($opt->names) {
    open(my $nh, '<', $opt->names) || die "Could not open name file: $!";
    while (! eof $nh) {
        my $line = <$nh>;
        if ($line =~ /(\d+)\t(.+)/) {
            $nameH{$1} = $2;
        }
    }
}
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
# This maps each cell line to its panel type.
my %lineType;
# This contains the name of the current group.
my $groupName = '';
my $drugs = 0;
my $groups = 0;
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    my @cols = split /,/, $line;
    # Compute the drug names. Note drug 2 might not exist.
    my $drug1 = $nameH{$cols[8]} // $cols[8];
    my $drug2 = ($cols[14] ? ($nameH{$cols[14]} // $cols[14]) : '');
    my $group = "$drug1,$drug2,$cols[28]";
    ($cols[8],$cols[14]) = ($drug1, $drug2);
    # Save the panel/line mapping.
    $lineType{$cols[28]} = $cols[27];
    # Is this the end of a group?
    if ($group ne $groupName) {
        # Yes. Set up for the next group.
        if (! $opt->doubles || $drugs == 2) {
            ProcessGroup($groupName, $lines, \@saved, $criterion) if $groupName;
        }
        $groupName = $group;
        $drugs = ($cols[14] ? 2 : 1);
        $lines = [];
    }
    push @$lines, \@cols;
}
# We don't process the residual-- the file comes with a trailer line.
# This hash contains counts for our summary report.
my %summary;
# This hash contains pairs for our triples report.
my %pairs;
# This hash contains counts for the cell-lines report.
my %lineCounts;
# This hash contains panels found in the cell-lines report.
my %panels;
# If this is grouped output, resort it.
if ($opt->grouped) {
    print STDERR "Sorting groups.\n";
    my @temp = sort { grp_cmp($a->[2], $b->[2]) } @saved;
    @saved = @temp;
}
print STDERR "Printing groups.\n";
# Print the saved groups.
for my $groupData (@saved) {
    my ($score, $group, $name) = @$groupData;
    # Count the group.
    my ($drug1, $drug2, $cellLine) = split /,/, $name;
    $summary{$drug1}{$drug2}++;
    $summary{$drug2}{$drug1}++;
    if ($summary{$drug1}{$drug2} >= $opt->tripsig) {
        $pairs{$drug1}{$drug2} = 1;
    }
    # Tally the cell line's participation for the drugs.
    my $panel = $lineType{$cellLine};
    $lineCounts{$drug1}{$panel}++;
    $lineCounts{$drug2}{$panel}++;
    $panels{$panel} = 1;
    # Write out the group.
    for my $line (@$group) {
        print join("\t", $score, @$line) . "\n";
    }
}
# Print the summary report if requested.
if ($opt->summary) {
    print STDERR "Calculating summary.\n";
    # Compute the total and the report for each drug.
    my %totals;
    my %xref;
    for my $drug1 (keys %summary) {
        my $countH = $summary{$drug1};
        my @others = sort { $countH->{$b} <=> $countH->{$a} } keys %$countH;
        my $total = 0;
        for my $drug2 (@others) {
            my $n = $countH->{$drug2};
            if ($n >= $opt->significant) {
                $xref{$drug1}{$drug2} = $n;
                $xref{$drug2}{$drug1} = $n;
            }
            $total += $n;
        }
        $totals{$drug1} = $total;
    }
    # Now output the summary report.
    print STDERR "Printing summary.\n";
    open(my $oh, '>', $opt->summary) || die "Could not open summary file: $!";
    my @drugs = sort { $totals{$b} <=> $totals{$a} } keys %xref;
    print $oh join("\t", "Drug", "Count", @drugs) . "\n";
    for my $drug1 (@drugs) {
        my @line = ($drug1, $totals{$drug1});
        for my $drug2 (@drugs) {
            push @line, $xref{$drug1}{$drug2} // '';
        }
        print $oh join("\t", @line) . "\n";
    }
}
# Print the triples report if required.
if ($opt->triples) {
    print STDERR "Printing triples.\n";
    open(my $oh, '>', $opt->triples) || die "Could not open triples file: $!";
    print $oh join("\t", 'A', 'B', 'C', 'A with B', 'A with C', 'B with C') . "\n";
    for my $drug1 (sort keys %pairs) {
        my @others = sort keys %{$pairs{$drug1}};
        for (my $i = 0; $i < scalar @others; $i++) {
            my $drug2 = $others[$i];
            for (my $j = $i + 1; $j < scalar @others; $j++) {
                my $drug3 = $others[$j];
                if ($pairs{$drug2}{$drug3}) {
                    print $oh join("\t", $drug1, $drug2, $drug3, $summary{$drug1}{$drug2}, $summary{$drug1}{$drug3}, $summary{$drug2}{$drug3}) . "\n";
                }
            }
        }
    }
}
# Print the cell-lines report if required.
if ($opt->panels) {
    print STDERR "Printing panel report.\n";
    open(my $oh, '>', $opt->panels) || die "Could not open panels file: $!";
    my @panels = sort keys %panels;
    print $oh join("\t", 'drug', @panels) . "\n";
    # Sort the drugs by synergy.
    my %drugCounts;
    for my $drug (keys %lineCounts) {
        my $panelsH = $lineCounts{$drug};
        my $total = 0;
        for my $panel (keys %$panelsH) {
            $total += $panelsH->{$panel};
        }
        $drugCounts{$drug} = $total;
    }
    my @drugs = sort { $drugCounts{$b} <=> $drugCounts{$a} } keys %drugCounts;
    for my $drug (@drugs) {
        my $panelsH = $lineCounts{$drug};
        print $oh join("\t", $drug, map { $panelsH->{$_} // 0 } @panels) . "\n";
    }
}


# Sort function for grouped mode.
sub grp_cmp {
    my ($a, $b) = @_;
    my ($da, $ea, $ca) = split /,/, $a;
    my ($db, $eb, $cb) = split /,/, $b;
    my $retVal = ($da <=> $db) || ($ea <=> $eb) || ($ca cmp $cb);
    return $retVal;
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
        # Determine the basic criterion for keeping a group.
        if ($opt->min) {
            # Here we are keeping everything above the minimum.
            if ($score >= $opt->min) {
                @$saved = sort { $b->[0] <=> $a->[0] } (@$saved, [$score, $lines, $groupName]);
            }
        } else {
            # Keeping the top N. Sort the new group into the saved set.
            @$saved = sort { $b->[0] <=> $a->[0] } (@$saved, [$score, $lines]);
            # Insure we are not too big.
            if (scalar @$saved > $opt->n) {
                pop @$saved;
            }
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
