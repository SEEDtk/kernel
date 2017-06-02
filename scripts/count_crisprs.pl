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
use Shrub;
use ScriptUtils;
use SeedUtils;

=head1 Count CRISPR Components

    count_crisprs.pl [ options ] inDir

This script processes the output from L<find_crisprs.pl> to produce a count of the CRISPR components in the analyzed
genomes.

=head2 Parameters

The positional parameter is the name of the directory containing the C<.json> files produced by L<find_crisprs.pl>.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir',
        Shrub::script_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the input directory.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Input directory $inDir invalid or missing.";
}
# Get the genome JSON files.
opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
my @gFiles = grep { $_ =~ /^\d+\.\d+\.json$/ } readdir $dh;
closedir $dh;
# Track the totals in here.
my %totals;
# Print the headings.
print join("\t", qw(ID Name Spacers Arrays Min Max)) . "\n";
# Loop through the files.
for my $gFile (@gFiles) {
    # Extract the genome ID.
    my ($genomeID) = $gFile =~ /(\d+\.\d+)/;
    # Get its name.
    my ($name) = $shrub->GetFlat('Genome', 'Genome(id) = ?', [$genomeID], 'name');
    $name //= '** name not found';
    # Load the file.
    my $crisprList = SeedUtils::read_encoded_object("$inDir/$gFile");
    my $arrays = scalar @$crisprList;
    my ($spacers, $min, $max) = (0, 0, 0);
    for my $crispr (@$crisprList) {
        my $spacerList = $crispr->[3];
        my $count = scalar @$spacerList;
        $spacers += $count;
        if (! $min || $min > $count) { $min = $count; }
        if (! $max || $max < $count) { $max = $count; }
    }
    print join("\t", $genomeID, $name, $spacers, $arrays, $min, $max) . "\n";
    $totals{spacers} += $spacers;
    $totals{arrays} += $arrays;
    if ($min && $max) {
        if (! $totals{min} || $totals{min} > $min) { $totals{min} = $min; }
        if (! $totals{max} || $totals{max} < $max) { $totals{max} = $max; }
    }
}
print join("\t", 'TOTALS', '', $totals{spacers}, $totals{arrays}, $totals{min}, $totals{max}) . "\n";
