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
use SeedUtils;
use File::Copy::Recursive;
use Stats;
use GenomeTypeObject;

=head1 Compute GTO Genome Quality

    gto_quality.pl [ options ] inDir

This script will take as input a directory of L<GenomeTypeObject> files and run the SciKit quality check on each one.
The output will be a three-column file of genomeID, genome name, and quality measure.

=head2 Parameters

The single positional parameter is the input directory name.

The command-line options the following.

=over 4

=item temp

The name of a temporary working directory. The default is the SEEDtk temporary directory.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir',
        ['temp|t=s', 'temporary working directory', { default => $FIG_Config::temp }],
        );
# Compute the input directory.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid input directory $inDir.";
}
# Create a statistics object.
my $stats = Stats->new();
# Get a list of all the GTOs.
opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
my @gtos = sort grep { $_ =~ /^\d+\.\d+\.gto$/} readdir $dh;
closedir $dh;
print STDERR scalar(@gtos) . " GTOs found in $inDir.\n";
# Get the control options.
my $missing = $opt->missing;
# Get the minimum quality threshold.
my $min = $opt->min;
# Get the temporary directory.
my $tempDir = $opt->temp . "/gtoq_$$";
File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
eval {
    # Create the results directory.
    my $resultDir = "$tempDir/results";
    File::Copy::Recursive::pathmk($resultDir) || die "Could not create $resultDir: $!";
    # Loop through the input genomes.
    for my $gtoFile (@gtos) {
        print STDERR "Processing $gtoFile.\n";
        my $gtoFileName = "$inDir/$gtoFile";
        my $gto = GenomeTypeObject->new_from_file($gtoFileName);
        if ($gto->{domain} ne 'Bacteria' && $gto->{domain} ne 'Archaea') {
            # Here it is not prokaryotic.
            $stats->Add(genomeNotProk => 1);
        } else {
            # We have the genome. Get its ID and name.
            my $genomeID = $gto->{id};
            my $name = $gto->{scientific_name};
            # Clear the result directory.
            File::Copy::Recursive::pathrmdir($resultDir) || die "Could not clear $resultDir: $!";
            # Compute the quality.
            my $cmd = "gto_consistency $gtoFileName $resultDir $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $FIG_Config::global/roles.to.use";
            SeedUtils::run($cmd);
            # Read in the results.
            open(my $qh, "<$resultDir/evaluate.log") || die "Could not open $genomeID quality log: $!";
            my $quality;
            while (! eof $qh && ! defined $quality) {
                my $line = <$qh>;
                if ($line =~ /Consistency=\s+(\d+(?:\.\d+)?)%/) {
                    $quality = $1;
                }
            }
            # Get a quality range for statistical purposes.
            my $qType = int($quality/10) . "X";
            $stats->Add("quality$qType" => 1);
            print join("\t", $genomeID, $name, $quality) . "\n";
        }
    }
};
if ($@) {
    print STDERR "FATAL ERROR: $@\n";
}
print STDERR "Cleaning up $tempDir.\n";
File::Copy::Recursive::pathrmdir($tempDir);
print STDERR "Statistics:\n" . $stats->Show();