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

=head1 Parse Drug Name File

    parse_drug_name_file.pl [ options ] outDir

This script parses the input files and outputs two tab-delimited
files: C<drug_name>, which maps drug numbers to drug names, and C<name_drug> which maps drug
names to NSC numbers. The drug numbers will be prefixed by the source name (e.g. C<NSC.>) in order to match the
PATRIC convention.

=head2 Parameters

The positional parameter is the name of the output directory into which the two output files should
be placed.

The command-line options are used to specify the input files. At least one must be specified.

=over 4

=item NSC

Name of the file containing the NSC drug IDs (which should be the chemnames file downloaded from
L<https://wiki.nci.nih.gov/display/NCIDTPdata/Chemical+Data>).

=item GDSC

Name of the file containing the GDSC drug IDs (which should be a tab-delimited version of thethe Screened_Compounds
spreadsheet from L<http://www.cancerrxgene.org/downloads>).

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        ['NSC|nsc=s', 'NSC drug mapping file'],
        ['GDSC|gdsc=s', 'GDSC drug mapping file']
        );
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Output directory $outDir is missing or invalid.";
}
# Create the stats object.
my $stats = Stats->new();
# Create a hash of the input files.
my %types;
if ($opt->nsc) {
    $types{NSC} = $opt->nsc;
}
if ($opt->gdsc) {
    $types{GDSC} = $opt->gdsc;
}
# Open the output files.
print "Creating output files.\n";
open(my $nh, ">$outDir/name_drug") || die "Could not open name_drug file: $!";
open(my $dh, ">$outDir/drug_name") || die "Could not open drug_name file: $!";
# Loop through the file types.
for my $type (sort keys %types) {
    # This will remember the drug numbers found so far. We choose the first name for any given number (rather arbitrarily)
    # to be the main name.
    my %found;
    # Open the input file.
    print "Accessing $type input file.\n";
    open(my $ih, '<', $types{$type}) || die "Could not open $type input: $!";
    # Compute the match pattern.
    my $pattern;
    if ($type eq 'NSC') {
        $pattern = qr/^(\d+)\|([^|]+)/;
    } elsif ($type eq 'GDSC') {
        $pattern = qr/^(\d+)\t([^\t]+)/;
    }
    # This will count the records read for this file.
    my $count = 0;
    # Loop through the file.
    while (! eof $ih) {
        my $line = <$ih>;
        $count++;
        $stats->Add(lineIn => 1);
        if ($line =~ $pattern) {
            $stats->Add(lineParsed => 1);
            # Parse out the drug number and name.
            my ($drug, $name) = ($1, $2);
            if (! $found{$drug}) {
                # Here it's a new drug. Output the number/name map.
                $stats->Add(drugsFound => 1);
                print $dh "$type.$drug\t$name\n";
                $found{$drug} = 1;
            }
            # Always output the name/number map.
            print $nh "$name\t$type.$drug\n";
            $stats->Add(namesFound => 1);
        }
        print "$count lines processed.\n" if $count % 5000 == 0;
    }
    print "$count lines in $type file.\n";
}

print "All done.\n" . $stats->Show();