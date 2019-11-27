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

=head1 Compare Genome to Close_Anno.pl Output

    p3x-genome-function-file.pl [options] genomeID

This is a special-purpose script that compares the output of L<Close_Anno.pl> to an existing genome.  It downloads all the
CDS features of the genome and analyzes the stop locations and roles.  The resulting tables are compared to the genome
information in the standard input.

=head2 Parameters

The positional parameter is the name of the genome ID.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The standard input is tab-delimited with headers.  Each record contains (0) a contig ID, (1) a start location, (2) a
strand, (3) a stop location, and (4) a functional assignment.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use RoleParse;
use SeedUtils;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('genomeID', P3Utils::ih_options(),
    );
my $stats = Stats->new();
# Get the genome ID.
my ($genomeID) = @ARGV;
if (! $genomeID) {
    die "No genome ID specified.";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# This hash will count the roles.  It gets +1 for each role in the genome, and -1 for each role in the input file.
# The result is the number of missing roles (negative for extra roles).  The key is role checksum.
my %roleCounts;
# This hash maps role checksums to role names.
my %checksum;
# This hash maps each stop location in the genome to a [length, roles] list.  A stop location is represented as
# contigID-dir-stop.
my %orfs;
# Loop through the features in the PATRIC genome.
print STDERR "Analyzing genome $genomeID.\n";
my $features = P3Utils::get_data($p3, feature => [['eq', 'genome_id', $genomeID], ['eq', 'feature_type', 'CDS']],
        ['sequence_id', 'start', 'end', 'strand', 'product']);
my ($fcount, $ftotal) = (0, scalar(@$features));
print STDERR "$ftotal features found.\n";
for my $feature (@$features) {
    $stats->Add(patricPeg => 1);
    my ($contig, $start, $stop, $dir, $function) = @$feature;
    # Compensate for strand.
    my $len = $stop - $start + 1;
    ($start,$stop) = ($stop,$start) if $dir eq '-';
    # Parse the function.
    my $roles = process_function($function, 1);
    # Store the ORF in the ORF hash.
    $orfs{"$contig$dir$stop"} = [$len, @$roles];
    $fcount++;
    print STDERR "$fcount of $ftotal PATRIC features processed.\n" if ($fcount % 100 == 0);
}
# Open the input file.
my $ih = P3Utils::ih($opt);
# Skip the header.
my $line = <$ih>;
# Loop through the input file.
print STDERR "Analyzing input file.\n";
$fcount = 0;
while (! eof $ih) {
    my ($contig, $start, $dir, $stop, $function) = P3Utils::get_fields($ih);
    $stats->Add(newPeg => 1);
    # Compute the ORF.
    my $orf = "$contig$dir$stop";
    # Compute the length.
    my $len = $stop - $start;
    $len = -$len if $len < 0;
    $len++;
    # Parse the function.
    my $roles = process_function($function, -1);
    if (! $orfs{$orf}) {
        $stats->Add(extraOrf => 1);
    } else {
        my ($len2, @roles2) = @{$orfs{$orf}};
        if ($len2 > $len) {
            $stats->Add(patricLonger => 1);
        } elsif ($len > $len2) {
            $stats->Add(newLonger => 1);
        } else {
            $stats->Add(sameLength => 1);
        }
        if (scalar @$roles > scalar @roles2) {
            $stats->Add(newMoreRoles => 1);
        } elsif (scalar @$roles < scalar @roles2) {
            $stats->Add(patricMoreRoles => 1);
        } else {
            for (my $i = 0; $i < @roles2; $i++) {
                if ($roles->[$i] eq $roles2[$i]) {
                    $stats->Add(roleMatch => 1);
                } else {
                    $stats->Add(roleDifferent => 1);
                }
            }
        }
    }
    $fcount++;
    print STDERR "$fcount new features processed.\n" if ($fcount % 100 == 0);
}
# Unspool the role comparisons.
print "role\tcount\tpatric\tnew\n";
my @roles = sort { $checksum{$a} cmp $checksum{$b} } keys %roleCounts;
for my $checksum (@roles) {
    my ($count, $patric, $new) = @{$roleCounts{$checksum}};
    if ($count) {
        print "$checksum{$checksum}\t$count\t$patric\t$new\n";
        if ($count > 0) {
            $stats->Add(patricExtraRole => $count);
        } else {
            $stats->Add(newExtraRole => -$count);
        }
    }
}
print STDERR "All done.\n" . $stats->Show();

## Process a function.  The function is parsed into roles and stored in the role hashes.  The type
## is +1 for PATRIC and -1 for NEW.  The return value is a reference to a list of the checksums.
sub process_function {
    my ($function, $type) = @_;
    my @retVal;
    my $kind = ($type == 1 ? 'patric' : 'new');
    my $idx = ($type == 1 ? 2 : 1);
    my @roles = SeedUtils::roles_of_function($function);
    for my $role (@roles) {
        $stats->Add($kind . 'Roles' => 1);
        my $checksum = RoleParse::Checksum($role);
        if (! $checksum{$checksum}) {
            $stats->Add(uniqueRoles => 1);
            $checksum{$checksum} = $role;
            if ($type == -1) {
                $stats->Add(newUniqueRole => 1);
            }
            $roleCounts{$checksum} = [0,0,0];
        }
        $roleCounts{$checksum}[0] += $type;
        $roleCounts{$checksum}[$idx]++;
        push @retVal, $checksum;
    }
    return [ sort @retVal ];
}