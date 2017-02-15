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
The output will be a three-column file of genomeID, genome name, coarse quality measure, and fine quality measure.

=head2 Parameters

The single positional parameter is the input directory name.

The command-line options the following.

=over 4

=item list

If specified, a list of GTOs to process. The file should be tab-delimited, with genome IDs in the first column.
The files processed will be those with a name equal to the genome ID and a suffix of C<.gto> (e.g. C<100226.1.gto>
for C<100226.1>), and all must be in the input directory.

=item temp

The name of a temporary working directory. The default is the SEEDtk temporary directory.

=item keep

If specified, the output logs will be kept. Do not do this for large runs.

=item filter

A filter file of roles to use. This is passed as a parameter to the L<gto_consistency.pl> script. If C<0> is
specified, then no filtering is used. The file should be tab-delimited, with role IDs in the first column. It
must be the same filter file used to generate the function predictors with L<build_matrix.pl>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir',
        ['list|f=s', 'file containing list of GTOs to process'],
        ['temp|t=s', 'temporary working directory', { default => $FIG_Config::temp }],
        ['filter=s', 'filter file of roles to use, or 0 for no filtering', { default => "$FIG_Config::global/roles.to.use"}],
        ['keep', 'keep output logs'],
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
my @gtos;
my $list = $opt->list;
if ($list) {
    open(my $ih, "<$list") || die "Could not open $list: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\d+\.\d+)/) {
            push @gtos, "$1.gto";
        }
    }
    print STDERR scalar(@gtos) . " GTOs found in file $list.\n";
} else {
    opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
    @gtos = sort grep { $_ =~ /^\d+\.\d+\.gto$/} readdir $dh;
    closedir $dh;
    print STDERR scalar(@gtos) . " GTOs found in directory $inDir.\n";
}
# Compute the filter file.
my $roles_to_use = $opt->filter;
if (! $roles_to_use) {
    $roles_to_use = '';
}
# Get the temporary directory.
my $tempDir = $opt->temp . "/gtoq_$$";
File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
eval {
    # Create the results directory.
    File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
    # Loop through the input genomes.
    for my $gtoFile (@gtos) {
        print STDERR "Processing $gtoFile.\n";
        my $gtoFileName = "$inDir/$gtoFile";
        my $gto = GenomeTypeObject->create_from_file($gtoFileName);
        if ($gto->{domain} ne 'Bacteria' && $gto->{domain} ne 'Archaea') {
            # Here it is not prokaryotic.
            $stats->Add(genomeNotProk => 1);
        } else {
            # We have the genome. Get its ID and name.
            my $genomeID = $gto->{id};
            my $name = $gto->{scientific_name};
            # Clear the result directory.
            my $resultDir = "$tempDir/results";
            # If we are keeping the logs, we need a separate directory for each GTO.
            if ($opt->keep) {
                $resultDir = "$tempDir/$gtoFileName";
            }
            File::Copy::Recursive::pathrmdir($resultDir) || die "Could not clear $resultDir: $!";
            # Compute the quality.
            my $cmd = "gto_consistency $gtoFileName $resultDir $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $roles_to_use";
            SeedUtils::run($cmd);
            # Read in the results.
            open(my $qh, "<$resultDir/evaluate.log") || die "Could not open $genomeID quality log: $!";
            my ($cquality, $fquality);
            while (! eof $qh) {
                my $line = <$qh>;
                if ($line =~ /Coarse_Consistency=\s+(\d+(?:\.\d+)?)%/) {
                    $cquality = $1;
                } elsif ($line =~ /Fine_Consistency=\s+(\d+(?:\.\d+)?)%/) {
                    $fquality = $1;
                }
            }
            # Get a quality range for statistical purposes.
            my $qType = int($fquality/10) . "X";
            $stats->Add("Fquality$qType" => 1);
            $qType = int($cquality/10) . "X";
            $stats->Add("Cquality$qType" => 1);
            print join("\t", $genomeID, $name, $cquality, $fquality) . "\n";
        }
    }
};
if ($@) {
    print STDERR "FATAL ERROR: $@\n";
}
if ($opt->keep) {
    print STDERR "Results kept in $tempDir.\n";
} else {
    print STDERR "Cleaning up $tempDir.\n";
    File::Copy::Recursive::pathrmdir($tempDir);
}
print STDERR "Statistics:\n" . $stats->Show();