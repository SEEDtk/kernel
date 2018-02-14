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
use Stats;

=head1 Merge Drug Data From Multiple Sources

    merge_drug_groups.pl [ options ] workDir

This is a rather arcane script that mixes data from several sources into a final report. It accepts as input a work directory containing the following
tab-delimited files. Some files come from the Drug Synergy Prediction Machine Learning Project (RICK) and some from the L<get_drug_groups.pl> script (MIKE).
The RICK files contain information about drugs and cell lines in an unusual format-- a drug is specified by an NSC drug ID (e.g. C<NSC.26980>) followed by
the drug name in parentheses. Spaces in the drug name are replaced by underscores. A cell line is specified by the prefix C<GDSC.> and then the cell
line name. When a drug pair is specified, the drug specifiers are separated by lower-case C<x>. If a drug pair is specified with a cell line, the cell
line is first and once again lower-case C<x> is used as a separator.

We will extract the data relating to NCI60 cell lines and specific drug pairs. The data from the original RICK detail file will be combined with growth
data and pair statistics from the MIKE files to permit easy comparison.

Note that not all columns of the input file are used. The index numbers below are array indices, so the first column is column 0. All files have header
lines.

=over 4

=item RICK files

=over 8

=item rick_pair_counts.tbl

Contains (0) the synergy rank (0-based), (1) the drug pairs of interest, (2) the synergy count, and (3) the mean synergy.

=item rick_growth_data.tbl

Contains (1) the drug pair and cell line, (2) the mean predicted growth, and (10) the expected growth.

=back

=item MIKE files

=over 8

=item pairs.tbl

Contains (0) the first drug name, (1) the second drug name, (2) the mean percent growth, (3) the number of significant synergies, and (4) the
mean percent growth for cell lines with significant synergy.

=item growth.csv

This is a comma-delimited file, downloaded from the NCI60 website.

Contains (8) the first drug number, (11) the first drug concentration, (14) the second drug number, (17) the second drug concentration, (19) the
observed percent growth, (24) the expected percent growth, and (28) the cell line name.

=back

=item Support files

=item names.tbl

Maps drug numbers to drug names. This file has no headers, and the columns are (0) NSC drug number, and (1) drug name.

=back

=head2 Parameters

The sole positional parameter should be the name of a directory containing the files described above.

=head2 Output Files

There are two output files.

=over 4

=item combo_pairs.tbl

For each drug pair, contains (0) the first drug name, (1) the second drug name, (2) the ranking in the MIKE files by mean percent growth,
(3) the ranking in the RICK files.

=item combo_detail.tbl

For each drug pair, cell line, and concentration combination in C<details.tbl>, contains (0) the first drug name, (1) the second drug name, (2) the
cell line name, (3) the RICK mean predicted growth, (4) the RICK expected growth, (5) the first drug concentration, (6) the second drug concentration,
(7) the observed percent growth, and (8) the expected percent growth.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir'
        );
my $stats = Stats->new();
# Get the working directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    die "Invalid or missing work directory $workDir.";
} elsif (! -s "$workDir/rick_pair_counts.tbl" || ! -s "$workDir/rick_growth_data.tbl") {
    die "Working directory $workDir is missing one of the RICK files.";
} elsif (! -s "$workDir/pairs.tbl" || ! -s "$workDir/details.tbl") {
    die "Working directory $workDir is missing one of the MIKE files.";
}
# Read in the name file.
open(my $ih, "<$workDir/names.tbl") || die "Could not open names.tbl: $!";
my %names;
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /(\d+)\t(.+)/) {
        $names{$1} = $2;
        $stats->Add(drugNameRead => 1);
    }
}
print scalar(keys %names) . " drug names read.\n";
close $ih; undef $ih;
# Open the file containing the pairs of interest.
open($ih, "<$workDir/rick_pair_counts.tbl") || die "Could not open rick_pair_counts.tbl: $!";
# Throw away the header line.
my $header = <$ih>;
# This is a hash keyed on drug1\tdrug2. The value is the RICK rank.
my %pairs;
print "Reading rick_pair_counts.tbl.\n";
# Loop through the pairs.
while (! eof $ih) {
    my ($rank, $pair) = get_cols($ih, 0, 1);
    my ($d1, $d2) = get_drugs($pair);
    if (! $d1) {
        $stats->Add(rPairBadLine => 1);
    } else {
        $stats->Add(rPairGoodLine => 1);
        $pairs{"$d1\t$d2"} = $rank + 1;
    }
}
print scalar(keys %pairs) . " drug pairs found.\n";
close $ih; undef $ih;
# Now we read the MIKE details table. This tells us which cell lines we care about.
# We'll keep them in here.
my %lines;
# This is a 2D hash that maps drug1\tdrug2 and cell line to a list containing [concentration1, concentration2, growth, expected-growth] tuples.
my %growths;
# Loop through the details file.
open($ih, "<$workDir/growth.csv") || die "Could not open growth.csv: $!";
print "Reading growth.csv\n";
# Throw away the header line.
$header = <$ih>;
my $count = 0;
# Loop through the details.
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my @fields = split /,/, $line;
    my ($drug1, $conc1, $drug2, $conc2, $pctGrowth, $ePctGrowth, $cline) = @fields[8, 11, 14, 17, 19, 24, 28];
    if (! $names{$drug1} || ! $names{$drug2}) {
        # This indicates a solo-drug situation. We skip it.
        $stats->Add(detailBadLine => 1);
    } else {
        my $drugPair = join("\t", sort ($names{$drug1}, $names{$drug2}));
        $count++; print "$count lines processed.\n" if $count % 10000 == 0;
        # Is this a pair of interest?
        if (! $pairs{$drugPair}) {
            $stats->Add(detailSkipLine => 1);
        } else {
            # Yes. Memorize the cell line.
            $lines{$cline}++;
            # Check for existing growth data.
            my $list = $growths{$drugPair}{$cline};
            if (! $list) {
                $list = [];
                $growths{$drugPair}{$cline} = $list;
                $stats->Add(detailTriplets => 1);
            }
            # Store the growth data.
            push @$list, [$conc1, $conc2, $pctGrowth, $ePctGrowth];
            $stats->Add(detailUseLine => 1);
        }
    }
}
close $ih; undef $ih;
print scalar(keys %lines) . " cell lines found.\n";
# Validate the drug pairs.
my @allPairs = sort keys %pairs;
for my $drugPair (@allPairs) {
    if (! $growths{$drugPair}) {
        print "No records in growth.csv found for $drugPair.\n";
        $stats->Add(missingDrugPair => 1);
        delete $pairs{$drugPair};
    }
}
# Now we need the RICK growth data.
open($ih, "<$workDir/rick_growth_data.tbl") || die "Could not open rick_growth_data.tbl: $!";
print "Reading rick_growth_data.tbl.\n";
# This is a 2D hash that maps drug pair and cell line to [predicted growth, expected growth].
my %rGrowths;
# Throw away the header line.
$header = <$ih>;
# Loop through the file.
my %clFound = map { $_ => 0 } keys %lines;
my %clAll = map { $_ => 0 } keys %lines;
$count = 0;
while (! eof $ih) {
    my ($label, $predGrowth, $expGrowth) = get_cols($ih, 1, 2, 10);
    my ($cl, $d1, $d2) = get_label($label);
    if (! $cl) {
        $stats->Add(growthBadLine => 1);
    } else {
        $clAll{$cl}++;
        my $drugPair = "$d1\t$d2";
        if (! $pairs{$drugPair} || ! $lines{$cl}) {
            $stats->Add(growthSkipLine => 1);
        } else {
            # Here we have a valid drug pair and cell line.
            $clFound{$cl}++;
            $rGrowths{$drugPair}{$cl} = [$predGrowth * 100, $expGrowth * 100];
            $count++;
            $stats->Add(growthUseLine => 1);
        }
    }
}
close $ih; undef $ih;
print "$count triples found in file.\n";
print "Cell line analysis\nLine\tRick\tMike\tRick.All\n";
for my $cl (sort keys %lines) {
    print "$cl\t$clFound{$cl}\t$lines{$cl}\t$clAll{$cl}\n";
}
print "\n";
# Finally, the MIKE pair data. We need to sort these to compute ranks.
print "Reading pairs.tbl.\n";
open($ih, "<$workDir/pairs.tbl") || die "Could not open pairs.tbl: $!";
my %ranks;
# Throw away the header line.
$header = <$ih>;
# Loop through the file.
$count = 0;
while (! eof $ih) {
    my ($d1, $d2, $pctGrowth) = get_cols($ih, 0, 1, 2);
    my $drugPair = join("\t", sort ($d1, $d2));
    if (! $pairs{$drugPair}) {
        $stats->Add(pairSkipLine => 1);
    } else {
        $ranks{$drugPair} = $pctGrowth;
        $stats->Add(pairUseLine => 1);
    }
}
close $ih; undef $ih;
print "Sorting pairs.\n";
my @pairs = sort { $ranks{$b} <=> $ranks{$a} } keys %ranks;
my $r = 1;
for my $pair (@pairs) {
    $ranks{$pair} = $r;
    $r++;
}
print "All data read in.\n";
# Write the combo pairs file.
print "Writing combo_pairs.tbl.\n";
open(my $oh, ">$workDir/combo_pairs.tbl") || die "Could not open combo_pairs.tbl: $!";
# Write the header.
print_cols($oh, qw(drug1 drug2 rank pred_rank));
# Loop through the pairs.
for my $pair (@pairs) {
    print_cols($oh, $pair, $ranks{$pair}, $pairs{$pair});
    $stats->Add(pairsOut => 1);
}
close $oh; undef $oh;
# Write the combo details file.
print "Writing combo_details.tbl.\n";
open($oh, ">$workDir/combo_details.tbl") || die "Could not open combo_details.tbl: $!";
# Write the header.
print_cols($oh, qw(drug1 drug2 cline mean_pred_growth mean_exp_growth conc1 conc2 growth exp_growth));
# Loop through the lines.
for my $pair (sort keys %growths) {
    my $lineHash = $growths{$pair};
    my $rickHash = $rGrowths{$pair};
    if (! $rickHash) {
        $stats->Add(unmatchedPair => 1);
    } else {
        $stats->Add(matchedPair => 1);
        for my $line (sort keys %$lineHash) {
            my $lineRickData = $rickHash->{$line};
            if (! $lineRickData) {
                $stats->Add(unmatchedLine => 1);
            } else {
                my ($meanPredGrowth, $meanExpGrowth) = @$lineRickData;
                my $concData = $lineHash->{$line};
                for my $concDatum (@$concData) {
                    print_cols($oh, $pair, $line, $meanPredGrowth, $meanExpGrowth, @$concDatum);
                }
            }
        }
    }
}
close $oh; undef $oh;
print "All done.\n" . $stats->Show();


# Print the specified columns.
sub print_cols {
    my ($oh, @cols) = @_;
    print $oh join("\t", @cols) . "\n";
}

# Get the listed columns from the next input line.
sub get_cols {
    my ($ih, @cols) = @_;
    my $line = <$ih>;
    chomp $line;
    my @fields = split /\t/, $line;
    my @retVal = map { $fields[$_] } @cols;
    return @retVal;
}

# Get the drug names from a drug pair.
sub get_drugs {
    my ($drugPair) = @_;
    my ($d1, $d2) = ($drugPair =~ /^NSC\.(\d+)\(.+\)xNSC\.(\d+)\(.+\)/);
    my @retVal;
    if ($d1) {
        $d1 = $names{$d1} // '?';
        $d2 = $names{$d2} // '?';
        @retVal = sort ($d1, $d2);
    }
    return @retVal;
}

# Get the cell line and drug names from an experiment label.
sub get_label {
    my ($label) = @_;
    my ($cl, $d1, $d2) = ($label =~ /^GDSC\.([^x]+)xNSC\.(\d+)\(.+\)xNSC\.(\d+)\(.+\)/);
    my @retVal;
    if ($cl) {
        $d1 = $names{$d1} // '?';
        $d2 = $names{$d2} // '?';
        @retVal = ($cl, sort ($d1, $d2));
    }
    return @retVal;
}
