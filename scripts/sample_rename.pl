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
use GPUtils;

=head1 Rename Samples in Genome Packages

    sample_rename.pl [ options ] pkgDir

This script reads an input file with sample name mappings and modifies the data files of the packages in the specified
package directory to rename the samples.

=head2 Parameters

The positional parameter is the name of the genome package directory whose packages are to be modified.

The command-line options are those found in L<ScriptUtils/ih_options> (to select the input file).

The input file should be tab-delimited, the first column being genome IDs or old sample names and the second being new sample names.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('pkgDir', ScriptUtils::ih_options(),
        );
# Validate the parameters.
my ($pkgDir) = @ARGV;
if (! $pkgDir) {
    die "No package directory specified.";
} elsif (! -d $pkgDir) {
    die "Invalid or missing package directory $pkgDir.";
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Read the mapping.
my %sampleMap;
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /(\S+)\t(\S+)/) {
        $sampleMap{$1} = $2;
        print "Mapping $1 to $2.\n";
    }
}
# Read the package directory.
print "Reading $pkgDir.\n";
my $gHash = GPUtils::get_all($pkgDir);
# Loop through the packages.
for my $genome (sort keys %$gHash) {
    print "Processing $genome.\n";
    my $dataHash = GPUtils::get_data($gHash, $genome);
    my $oldName = $dataHash->{'Sample Name'};
    my $newName = $sampleMap{$genome} // $sampleMap{$oldName};
    if ($newName) {
        $dataHash->{'Sample Name'} = $newName;
        print "Changing $oldName to $newName.\n";
        # Rewrite the data file.
        my $genomeDir = $gHash->{$genome};
        open(my $oh, '>', "$genomeDir/data.tbl") || die "Could not open $genome data table file: $!";
        for my $key (sort keys %$dataHash) {
            print $oh "$key\t$dataHash->{$key}\n";
        }
    }
}
## TODO process the input to produce the output.