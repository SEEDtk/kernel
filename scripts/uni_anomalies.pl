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

=head1 Display Universal Role Anomalies.

    uni_anomalies.pl [ options ] uniFile

This script runs through a list of universal roles. For each, it lists the well-behaved genomes
that do not contain the role. For each genome, the name, ID, and dna-size are listed.

=head2 Parameters

The positional parameter is the name of a tab-delimited file containing universal role IDs in the first column.
If the parameter is omitted, the list of universal roles will be read from the database.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('uniFile',
                Shrub::script_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the universal role list.
my %unis;
my ($uniFile) = @ARGV;
if (! $uniFile) {
    %unis = map { $_ => 1 } $shrub->GetFlat('Function', 'Function(universal) = ?', [1], 'id');
} elsif (! -s $uniFile) {
    die "$uniFile missing or empty.";
} elsif (! open(my $uh, "<$uniFile")) {
    die "Could not open $uniFile: $!";
} else {
    while (! eof $uh) {
        my $line = <$uh>;
        if ($line =~ /^(\S+)/) {
            $unis{$1} = 1;
        }
    }
}
# Get all the genomes.
my %gHash = map { $_->[0] => [$_->[1], $_->[2]] } $shrub->GetAll('Genome', 'Genome(well-behaved) = ?',
        [1], 'id name dna-size');
# Get a sorted list of genome IDs.
my @genomes = sort keys %gHash;
# Compute the number of expected genomes containing functions.
print scalar(@genomes) . " well-annotated genomes found.\n";
# Count the number of times a genome was not found.
my %unFound;
# Count the universal proteins.
my $uniCount = 0;
# Loop through the universal functions.
for my $role (sort keys %unis) {
    my ($roleDesc) = $shrub->GetFlat('Function', 'Function(id) = ?', [$role], 'description');
    if (! $roleDesc) {
        die "$role not found in database.";
    }
    print "\nGenomes missing $role: $roleDesc.\n";
    $uniCount++;
    # Get all the genomes containing the functions.
    my %found = map { $_ => 1 } $shrub->GetFlat('Role2Function Function2Feature Feature2Genome',
            'Role2Function(from-link) = ? AND Function2Feature(security) = ?',
            [$role, 2], 'Feature2Genome(to-link)');
    # Loop through the full list of useful genomes, looking for any that were not found above.
    for my $genome (@genomes) {
        if (! $found{$genome}) {
            print join("\t", '', $genome, @{$gHash{$genome}}) . "\n";
            $unFound{$genome}++;
        }
    }
}
my $threshhold = $uniCount/2;
my @unFounded = sort { $unFound{$b} <=> $unFound{$a} } keys %unFound;
print "\n\nSuspicious genomes (threshhold $threshhold).\n";
for my $genome (@unFounded) {
    print join("\t", '', $unFound{$genome}, $genome, @{$gHash{$genome}}) . "\n";
}