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
use File::Copy::Recursive;

=head1 Produce Cleaned FASTA Files for Packages in a Directory

    clean_packages.pl [ options ] packageDir

This script runs through the packages in a directory producing improved FASTA files that can then be fed to L<package_cleaned_bins.pl> to produce
(hopefully) improved binning results. For each package, it will perform the following steps.

=over 4

=item 1

Run L<gto_consistency.pl> using the C<bin.gto> file and the C<FunctionPredictors.Big> predictor set to compute the good and bad roles. The output
will be put into a temporary directory.

=item 2

Run L<analyze_consistency.pl> to compute the good and bad contigs. The output will be put into the package directory.

=item 3

Run L<improve_consistency.pl> to remove the bad contigs. This will create a new C<bin>I<X>C<.fa> file in the package directory, where I<X> is the
specified suffix.

=back

After processing, L<package_cleaned_bins.pl> can be used to create new packages from the new FASTA files.

=head2 Parameters

The positional parameter is the name of the directory containing the packages.

The command-line options are the following.

=over 4

=item missing

If specified, only packages that have not been previously cleaned using the specified suffix will be processed. In other words, in the default case
(suffix is C<1>), packages with a C<bin1.fa> will not be cleaned.

=item force

If specified, then all packages will be cleaned. Otherwise, packages already marked as good will be skipped.

=item suffix

The suffix to use when creating the new FASTA files. The default is C<1>, which will create C<bin1.fa> files.

=item predictor

Name of the function predictor directory to use. The default is C<FunctionPredictors.Big> in the Data directory.

=item temp

Name of the temporary directory to use for evaluation output. The default is generated in the Temp directory.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir',
        ['missing', 'only process new packages'],
        ['suffix=s', 'suffix to use when creating new FASTA files', { default => '1' }],
        ['force', 'process packages even if they are already good'],
        ['predictor=s', 'function predictor directory to use', { default => "$FIG_Config::data/FunctionPredictors.Big" }],
        ['temp=s', 'temporary output directory'],
        );
# Get the options.
my $missing = $opt->missing;
my $suffix = $opt->suffix;
my $predictor = $opt->predictor;
my $force = $opt->force;
# Check the predictor directory for a roles.to.use.
my $rolesToUse;
if (-s "$predictor/roles.to.use") {
    $rolesToUse = "$predictor/roles.to.use";
} else {
    $rolesToUse = "$FIG_Config::p3data/roles.to.use";
}
# Check the input directory.
my ($packageDir) = @ARGV;
if (! $packageDir) {
    die "No input directory specified.";
} elsif (! -d $packageDir) {
    die "Input directory $packageDir not found or invalid.";
}
# Create the temporary directory.
my $tempCreated;
my $temp = $opt->temp;
if (! $temp) {
    $temp = "$FIG_Config::temp/Eval$$";
    print "Creating temp directory $temp.\n";
    File::Copy::Recursive::pathmk($temp) || die "Could not create directory $temp: $!";
    $tempCreated = 1;
} elsif (! -d $temp) {
    print "Creating user temp directory $temp.\n";
    File::Copy::Recursive::pathmk($temp) || die "Could not create directory $temp: $!";
    mkdir $temp;
    $tempCreated = 1;
}
eval {
    # Get all the packages.
    print "Scanning $packageDir.\n";
    opendir(my $dh, $packageDir) || die "Could not open $packageDir: $!";
    my @packages = grep { -s "$packageDir/$_/bin.gto" } readdir $dh;
    close $dh; undef $dh;
    my $total = scalar @packages;
    print "$total packages found in $packageDir.\n";
    # Loop through the packages, processing each one.
    my $count = 0;
    for my $package (@packages) {
        $count++;
        print "Processing $package ($count of $total).\n";
        my $thisDir = "$packageDir/$package";
        my $skip;
        if (! $force && -s "$thisDir/quality.tbl") {
            open(my $ih, "$thisDir/quality.tbl") || die "Could not open quality.tbl: $!";
            my @fields = ScriptUtils::get_line($ih);
            if ($fields[16]) {
                print "Package marked as good-- skipping.\n";
                $skip = 1;
            }
        }
        if (! $skip && $missing && -s "$thisDir/bin$suffix.fa") {
            print "Cleaned FASTA already processed-- skipping.\n";
            $skip = 1;
        }
        if (! $skip) {
            my $rc = system("gto_consistency", "$thisDir/bin.gto", $temp, $predictor, "$FIG_Config::p3data/roles.in.subsystems", $rolesToUse);
            if ($rc) {
                die "Error in gto_consistency: rc = $rc.";
            }
            $rc = system("analyze_consistency", "$thisDir/bin.gto", $temp, $thisDir);
            if ($rc) {
                die "Error in analyze_consistency: rc = $rc.";
            }
            $rc = system("improve_consistency", "--input=$thisDir/contigs.tbl", "--suffix=$suffix", $thisDir);
            if ($rc) {
                die "Error in improve_consistency: rc = $rc.";
            }
            # Clear the temp directory.
            print "Clearing $temp.\n";
            File::Copy::Recursive::pathempty($temp) || die "Could not clear temp directory: $!";
        }
    }
};
# Save the error state;
my $error = $@;
# Clean up the temp directory.
if ($tempCreated) {
    print "Deleting temp directory $temp.\n";
    File::Copy::Recursive::pathrmdir($temp);
}
# Re-throw any error.
if ($error) {
    die "ERROR in processing: $error";
}