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
use GPUtils;
use Math::Round;

=head1 Find Significant Roles in Quality Evaluations

    significant_roles.pl [ options ] pDir1 pDir2 ... pDirN

This script analyzes all the quality evaluations for a set of genome packages and outputs the roles that made a difference.

=head2 Parameters

The positional parameters are the named of genome package directories to check.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('pDir1 pDir2 ... pDirN',
        );
my $stats = Stats->new();
# Get the package directories.
my (@pDirs) = @ARGV;
# This will track the number of times each role makes a difference.
my %rolesUsed;
# This will track the total genomes processed.
my $totalGenomes = 0;
for my $packageDir (@pDirs) {
    # Get the packages and loop through them.
    print STDERR "Processing $packageDir.\n";
    my $gHash = GPUtils::get_all($packageDir);
    for my $genome (keys %$gHash) {
        my $gDir = $gHash->{$genome};
        $stats->Add(packagesFound => 1);
        if (-s "$gDir/EvalBySciKit/evaluate.out") {
            print STDERR "Checking $genome.\n";
            $stats->Add(evaluationsFound => 1);
            open(my $ih, "<$gDir/EvalBySciKit/evaluate.out") || die "Could not open evaluation for $genome: $!";
            while (! eof $ih) {
                my ($id, $predicted, $actual) = ScriptUtils::get_line($ih);
                if (! defined $predicted || ! defined $actual) {
                    $stats->Add(badLine => 1);
                } else {
                    if ($predicted != $actual) {
                        $rolesUsed{$id}++;
                        $stats->Add(roleMattered => 1);
                    } else {
                        $stats->Add(roleNotMattered => 1);
                    }
                }
            }
            $totalGenomes++;
        }
    }
}
my @roleSort = sort { $rolesUsed{$b} <=> $rolesUsed{$a} } keys %rolesUsed;
print STDERR scalar(@roleSort) . " significant roles found.\n";
print join("\t", "ID\tCount\tPercent") . "\n";
for my $role (@roleSort) {
    my $uses = $rolesUsed{$role};
    my $pct = Math::Round::nearest(0.01, $uses * 100 / $totalGenomes);
    print join("\t", $role, $uses, $pct) . "\n";
}
print STDERR "All done.\n" . $stats->Show();

