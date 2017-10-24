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
use Stats;

=head1 Output Feature Coupling Profile Table

    get_feature_coupling_table.pl [ options ]

This script will output a file containing feature IDs, locations, roles, and (possibly) protein family IDs for a given list of genomes (found
in the input file).

=head2 Parameters

The input file should contain a list of genome IDs in the first column.

The command-line options are those found in L<Shrub/script_options> and L<ScriptUtils/ih_options> plus the following.

=over 4

=item header

Include headers in the output.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', Shrub::script_options(), ScriptUtils::ih_options(),
        ['headers|H', 'include a header line on the ouput']
        );
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Process the headers.
if ($opt->headers) {
    print join("\t", qw(feature contig location role family)) . "\n";
}
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
while (! eof $ih) {
    # Read the next genome ID.
    my $line = <$ih>;
    $stats->Add(linesIn => 1);
    if ($line =~ /(\d+\.\d+)/) {
        my $genomeID = $1;
        $stats->Add(genomes => 1);
        print STDERR "Processing $genomeID.\n";
        # We need to get the roles, the families, and the location data. The location data is first.
        print STDERR "Reading locations.\n";
        my %locs;
        my $q = $shrub->Get('Feature2Contig', 'Feature2Contig(from-link) LIKE ?', ["fig|$genomeID.peg.%"],
            'Feature2Contig(from-link) Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(len)');
        while (my $record = $q->Fetch()) {
            my ($fid, $contig, $start, $len) = $record->Values('Feature2Contig(from-link) Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(len)');
            $stats->Add(locRecord => 1);
            my $end = $start + $len - 1;
            my $oldLoc = $locs{$fid};
            # We merge the location data so that for each feature we have its overall start and end in a contig.
            if (! $oldLoc) {
                $locs{$fid} = [$contig, $start, $end];
            } else {
                my ($oldContig, $oldStart, $oldEnd) = @$oldLoc;
                if ($contig ne $oldContig) {
                    $stats->Add(multiContigFeature => 1);
                } else {
                    $oldLoc->[1] = $start if ($start < $oldStart);
                    $oldLoc->[2] = $end if ($end > $oldEnd);
                }
            }
        }
        # Now we get the roles for each feature. We only want single-role features.
        print STDERR "Processing roles.\n";
        my %roles;
        $q = $shrub->Get('Genome2Feature Feature2Function', 'Genome2Feature(from-link) = ? AND Feature2Function(security) = ?', [$genomeID, 0],
            'Feature2Function(from-link) Feature2Function(to-link)');
        while (my $record = $q->Fetch()) {
            $stats->Add(roleRecord => 1);
            my ($fid, $role) = $record->Values('Feature2Function(from-link) Feature2Function(to-link)');
            # Only process the role if it has no separators, that is, if it is a single-role function.
            if ($role =~ /[^a-z0-9]/i) {
                $stats->Add(multiRoleFidSkipped => 1);
            } else {
                push @{$roles{$fid}}, $role;
            }
        }
        # Finally, the families.
        print STDERR "Processing families.\n";
        my %families;
        $q = $shrub->Get('Genome2Feature Feature2Protein Protein2Family', 'Genome2Feature(from-link) = ?', [$genomeID],
            'Feature2Protein(from-link) Protein2Family(to-link)');
        while (my $record = $q->Fetch()) {
            $stats->Add(familyRecord => 1);
            my ($fid, $family) = $record->Values('Feature2Protein(from-link) Protein2Family(to-link)');
            push @{$families{$fid}}, $family;
        }
        # Now assemble the output.
        print STDERR "Assembling output.\n";
        for my $fid (sort keys %locs) {
            # Insure we have roles.
            my $roles = $roles{$fid};
            if (! $roles) {
                $stats->Add(fidNoRoles => 1);
                $roles = [''];
            }
            $stats->Add(fidOut => 1);
            # Format the feature and location data.
            my $loc = $locs{$fid};
            my ($contig, $start, $end) = @$loc;
            my @row = ($fid, $contig, "$start..$end");
            # Get any families there might be. If there are none, put in a blank.
            my $fams = $families{$fid};
            if (! $fams) {
                $fams = [''];
                $stats->Add(fidNoFamilies => 1);
            }
            if (scalar @$fams > 1) {
                $stats->Add(fidMultiFamily => 1);
            }
            for my $role (@$roles) {
                for my $fam (@$fams) {
                    print join("\t", @row, $role, $fam) . "\n";
                    $stats->Add(lineOut => 1);
                }
            }
        }
    }
}
print STDERR "All done.\n" . $stats->Show();
