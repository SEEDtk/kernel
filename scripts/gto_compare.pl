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
use GenomeTypeObject;
use Stats;

=head1 Compare Two GTOs

    gto_compare.pl [ options ] gto1 gto2

This script produces a report comparing the role profile of two L<GenomeTypeObject> instances. The GTOs must be
provided as files in JSON format. Each role will be converted to a Shrub role ID and counted. The counts are then
compared. A catchall will be used for roles that don't map to a shrub ID. Finally, there will be statistics on
the number of features and the DNA sequence length.

=head2 Parameters

The positional parameters are the names of the GTO files.  There must be at least two. All files must encode the
GTO in JSON format.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('gto1 gto2 ... gtoN',
                Shrub::script_options(),
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Create the statistics object.
my $stats = Stats->new();
# Get the GTO files.
my @gtoFiles = @ARGV;
my $gtoCount = scalar(@gtoFiles);
if ($gtoCount < 2) {
    die "$gtoCount arguments found. There must be at least two.";
}
# This hash counts role IDs.
my %roleCounts;
# This hash tracks feature IDs for each role.
my %roleFeats;
# This hash remembers the ID for each function.
my %funHash;
# This array counts the DNA.
my @dna;
# This array counts the features.
my @feats;
# This is the index of the current GTO.
my $i = 0;
# Loop through the GTO files.
for my $gtoFile (@gtoFiles) {
    # Get the actual GTO.
    if (! -s $gtoFile) {
        die "$gtoFile not found or empty.";
    }
    my $gto = GenomeTypeObject->create_from_file($gtoFile);
    print STDERR "Processing contigs of $gtoFile.\n";
    # Read through the contigs to get the DNA length.
    my $dnaLen = 0;
    my $contigsL = $gto->{contigs};
    for my $contig (@$contigsL) {
        my $len = length($contig->{dna});
        $dnaLen += $len;
        $stats->Add(contigs => 1);
        $stats->Add(dna => $len);
    }
    $dna[$i] = $dnaLen;
    # Read through the features to get the roles.
    print STDERR "Processing features of $gtoFile.\n";
    my $featCount = 0;
    my $featuresL = $gto->{features};
    for my $feature (@$featuresL) {
        $featCount++;
        $stats->Add(features => 1);
        # Compute the function ID.
        my $funID;
        my $function = $feature->{function};
        if (! $function) {
            $stats->Add(functionMissing => 1);
        } else {
            $stats->Add(functionRead => 1);
            if (exists $funHash{$function}) {
                # Here we already know the function.
                $funID = $funHash{$function};
                $stats->Add(functionFound => 1);
            } else {
                # Here we must get it from the database.
                $stats->Add(functionNotFound => 1);
                $funID = $shrub->desc_to_function($function);
                $funHash{$function} = $funID;
                if ($funID) {
                    $stats->Add(functionMapped => 1);
                } else {
                    $stats->Add(functionNotMapped => 1);
                }
            }
            if (! $funID) {
                $stats->Add(functionInvalid => 1);
                Increment(\%roleCounts, '(unknown)', $i);
            } else {
                my @roles = Shrub::roles_of_func($funID);
                for my $role (@roles) {
                    $stats->Add(roleProcessed => 1);
                    Increment(\%roleCounts, $role, $i);
                    push @{$roleFeats{$role}}, $feature->{id};
                }
            }
        }
    }
    # Record the feature count.
    $feats[$i] = $featCount;
    # Position forward in the GTO list.
    $i++;
}
# Loop through the role table, producing output.
for my $role (sort keys %roleCounts) {
    my $array = $roleCounts{$role};
    push @$array, 0 while (scalar(@$array) < $i);
    my $count = $array->[0];
    my $j = 1;
    while ($j < $i && $array->[$j] == $count) {
        $j++;
    }
    if ($j < $i) {
        # Here we want to print this role. Get its feature list.
        my $features = $roleFeats{$role};
        # If features exist, format the list.
        my $flist = ($features ? join(", ", sort @$features) : '');
        # Print the counts and the features.
        print join("\t", $role, @$array, $flist) . "\n";
        $stats->Add(roleMismatch => 1);
    } else {
        $stats->Add(roleMatch => 1);
    }
}
print "\n\n";
# Write the feature and DNA statistics.
print join("\t", '* Features', @feats) . "\n";
print join("\t", '* DNA', @dna) . "\n";
# Write the run statistics.
print STDERR "All done.\n" . $stats->Show();

## Increment an entry in an array inside a hash. Fill zeroes into the missing spaces.
sub Increment {
    my ($hash, $key, $i) = @_;
    # Get the array for the specified key.
    my $array;
    if (exists $hash->{$key}) {
        $array = $hash->{$key};
    } else {
        $array = [];
        $hash->{$key} = $array;
    }
    # Grow it until it is big enough to hold us.
    push @$array, 0 while (scalar(@$array) <= $i);
    # Increment the indexed item.
    $array->[$i]++;
}