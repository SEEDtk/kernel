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

=head1 Determine Roles to Use: Function Predictors Step 4

    build_roles_to_use.pl [ options ] probDir outDir

This script builds the C<roles.to.use> file in the specified probdir (which has been populated by L<build_predictor_set.pl>). This file is
then fed into L<build_matrix.pl> and the process is iterated again to produce a more stable prediction mechanism.

=head2 Parameters

The first positional parameter is the name of the probDir directory populated by the building of the predictor set. The second is the output
directory (which will be created if it does not yet exist). The C<roles.to.use> file will be built in the specified output directory.

The command-line options are as follows:

=over 4

=item min

Minimum acceptable accuracy for a predictor, in percent. The default is C<90.0>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('probDir outDir',
        ['min=f', 'minimum acceptable accuracy', { default => 90.0 }],
        );
my ($probDir, $outDir) = @ARGV;
if (! $probDir) {
    die "No probdir specified.";
} elsif (! -d $probDir) {
    die "Probdir $probDir missing or invalid.";
} elsif (! -d "$probDir/Predictors") {
    die "$probDir does not have a Predictors directory.";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
}
# Get the output file.
open(my $oh, ">$outDir/roles.to.use") || die "Could not open output file: $!";
# Loop through the predictors.
my $predDir = "$probDir/Predictors";
opendir(my $dh, $predDir) || die "Could not open Predictors: $!";
my @roles = sort grep { -f "$predDir/$_/Classifiers/RandomForestClassifier/accuracy" } readdir $dh;
closedir $dh;
my $total = scalar @roles;
print "$total predictors found.\n";
# Count the number of roles processed, kept, and rejected.
my ($count, $kept, $rejected) = (0,0,0);
# Loop through the predictors.
for my $role (@roles) {
    open(my $ih, "<$predDir/$role/Classifiers/RandomForestClassifier/accuracy") || die "Could not open accuracy file for $role: $!";
    my $line = <$ih>;
    chomp $line;
    my @values = split /\t/, $line;
    my (undef, $min, @data) = @values;
    for my $value (@data) {
        $min = $value if $value < $min;
    }
    if ($min < $opt->min) {
        print "$role rejected: min accuracy $min.\n";
        $rejected++;
    } else {
        print $oh join("\t", $role, $min) . "\n";
        $kept++;
    }
    $count++;
    print "$count of $total roles processed: $kept kept, $rejected rejected.\n";
}
# All done.
close $oh;
print "ALL DONE: $kept of $total kept, $rejected rejected.\n";
