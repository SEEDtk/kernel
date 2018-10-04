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
use File::Copy::Recursive;
use BinningReports;

=head1 Organize Role Deficiency Reports for Binning Result Packages

    organize_packages.pl [ options ] packageDir outDir

This script reads all packages in a package directory and organizes them by originating sample. For each sample,
an output file is produced containing a list of the excess roles and the missing roles based on the EvalG and EvalCon
results.

The output directory will have one sub-directory per sample, and each bin in the sample will have a report produced
named I<bin>C<.contigs.tbl>.  Each row of the report will contain a contig ID, the number of good roles in the
contig, and a list of the excess roles in the contig. The final row of the report will contain a hyphen (C<->) for
the contig ID, a count of C<0>, and a list of the missing roles in the bin.

=head2 Parameters

The positional parameters are the names of the input package directory and the name of an output directory to contain
the reports.

The command-line options are as follows.

=over 4

=item missing

If specified, reports will only be produced on samples for which a report does not already exist.

=item roleFile

The name of the file containing the master list of roles. This is a tab-delimited file with no headers. The first column of each
record is the role ID, the second is the role checksum, and the third is the role name. The default file is C<roles.in.subsystems>
in the global data directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir outDir',
        ['missing', 'only process packages for new samples'],
        ['roleFile|rolefile|roles=s', 'role file', { default => "$FIG_Config::global/roles.in.subsystems" }],
        );
# This hash will contain the samples already processed if MISSING is turned on.
my %doneSamples;
# Get the directories.
my ($packageDir, $outDir) = @ARGV;
if (! $packageDir) {
    die "No input directory specified.";
} elsif (! -d $packageDir) {
    die "Input directory $packageDir missing or invalid.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (-f $outDir) {
    die "Output directory $outDir invalid.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
}
# Create the role hash.
open(my $rh, '<', $opt->rolefile) || die "Could not open role file: $!";
# Loop through the roles.
my %roleMap;
while (! eof $rh) {
    my ($role, $checksum, $name) = ScriptUtils::get_line($rh);
    $roleMap{$checksum} = $role;
}
close $rh; undef $rh;
# Get all the bins.
print "Reading packages from $packageDir.\n";
my $gHash = GPUtils::get_all($packageDir);
# Loop through the bins.
for my $genome (sort keys %$gHash) {
    print "Processing $genome.\n";
    # Compute the sample.
    my $dataHash = GPUtils::get_data($gHash, $genome);
    my $sample = $dataHash->{'Sample Name'};
    my $sampDir = "$outDir/$sample";
    if (! -d $sampDir) {
        print "Creating $sampDir.\n";
        File::Copy::Recursive::pathmk($sampDir);
    }
    # Find out if we have already done this one.
    if ($opt->missing && -s "$sampDir/$genome.contigs.tbl") {
        print "Skipping-- $genome already processed.\n";
    } else {
        # Here it is worth analyzing the genome.
        my $gto = GPUtils::gto_of($gHash, $genome);
        BinningReports::UpdateGTO($gto, "$packageDir/EvalBySciKit", "$packageDir/EvalByCheckG", \%roleMap);
        # Now we need a directory of missing and excess roles.
        my (%missing, %excess);
        my $ppr = $gto->{genome_quality_measure}{problematic_roles_report};

    }
}
## TODO process the input to produce the output.