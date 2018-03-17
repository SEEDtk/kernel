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
use GenomeTypeObject;
use Stats;

=head1 Analyze Output from the Consistency Checker

    analyze_consistency.pl [ options ] gto evalDir

This script processes the output from the consistency checker for a specified L<GenomeTypeObject>. The output lists the problematic roles and
the contigs containing those roles. For each problematic role, we will output the role name, the predicted and actual occurrences, and the
features containing the role. For each contig, we will output the contig length, the number of features with good roles, and the features
containing problematic roles.

=head2 Parameters

The positional parameters are the file name of the GTO and the name of the directory into which the consistency-checker output was placed.

The output files will be placed into the specified directory:  C<contigs.tbl> will contain the contig analysis and C<roles.tbl> the
role analysis.

The command-line options will be as follows:

=over 4

=item roleFile

Name of the file containing the master ID-to-name mapping file for roles. The default is C<roles.in.subsystems> in the global data
directory. The file contains role IDs in the first column, checksums in the second column, and names in the last column.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('gto evalDir',
        ['roleFile=s', 'master role mapping file', { default => "$FIG_Config::global/roles.in.subsystems"}]
        );
# Create a statistics object.
my $stats = Stats->new();
# Check the positional parameters.
my ($gtoF, $evalDir) = @ARGV;
if (! $gtoF) {
    die "No GTO specified.";
} elsif (! -s $gtoF) {
    die "Missing or empty GTO file $gtoF.";
} elsif (! $evalDir) {
    die "Evaluation directory not specified.";
} elsif (! -s "$evalDir/evaluate.out") {
    die "$evalDir does not appear to contain evaluation output.";
}
# Read in the GTO.
print "Reading GTO from $gtoF.\n";
my $gto = GenomeTypeObject->create_from_file($gtoF);
# The first step is to read in the master list of roles. We need to map role IDs to role names.
my %roleNames;
print "Reading role names.\n";
open(my $ih, '<', $opt->rolefile) || die "Could not open roles.in.subsystems file: $!";
while (! eof $ih) {
    my ($id, $check, $name) = ScriptUtils::get_line($ih);
    $stats->Add(roleNameIn => 1);
    $roleNames{$id} = $name;
}
print scalar(keys %roleNames) . " role names found.\n";
close $ih; undef $ih;
# Now we need a map of role IDs to features and vice versa. This information can be found in roles.mapped.
my %roleFids;
my %fidRoles;
print "Reading role features.\n";
open($ih, "<$evalDir/roles.mapped") || die "Could not open roles.mapped: $!";
while (! eof $ih) {
    my ($name, $id, $fid) = ScriptUtils::get_line($ih);
    $stats->Add(roleMapIn => 1);
    push @{$roleFids{$id}}, $fid;
    push @{$fidRoles{$fid}}, $id;
}
print scalar(keys %roleFids) . " mapped roles found.\n";
print scalar(keys %fidRoles) . " features found.\n";
close $ih; undef $ih;
# Now we need to find the problematic roles. This hash maps each one to an [expected,actual] pair.
my %ppr;
# This hash remembers the good roles.
my %good;
print "Reading problematic roles.\n";
open($ih, "<$evalDir/evaluate.out") || die "Could not open evaluate.out: $!";
while (! eof $ih) {
    my ($id, $expect, $actual) = ScriptUtils::get_line($ih);
    $stats->Add(roleStatsIn => 1);
    if ($expect == $actual) {
        $good{$id} = 1;
        $stats->Add(roleGood => 1);
    } elsif ($expect != $actual) {
        $ppr{$id} = [$expect, $actual];
        $stats->Add(roleBad => 1);
    }
}
close $ih; undef $ih;
# Now we list the problematic roles.
open(my $oh, ">$evalDir/roles.tbl") || die "Could not open roles.tbl: $!";
print "Writing problematic roles.\n";
print $oh join("\t", qw(Role Description Expect Actual Fids)) . "\n";
for my $id (sort keys %ppr) {
    my $fids = $roleFids{$id} // [];
    my $tuple = $ppr{$id};
    print $oh join("\t", $id, $roleNames{$id}, @$tuple, @$fids) . "\n";
}
close $oh; undef $oh;
# Now we must analyze the contigs. This hash contains contig lengths.
my %contigLens;
# This hash contains the number of good roles.
my %contigGood;
# This hash contains a maps each contig to a list of the bad features.
my %contigBad;
# Get the list of contigs. First we compute the lengths and initialize the other hashes.
print "Reading contig list.\n";
my $contigList = $gto->{contigs};
for my $contig (@$contigList) {
    my $id = $contig->{id};
    my $len = length($contig->{dna});
    $contigLens{$id} = $len;
    $contigGood{$id} = 0;
    $contigBad{$id} = [];
    $stats->Add(contigIn => 1);
}
print scalar(keys %contigLens) . " contigs processed.\n";
# Now we need to loop through the features to place the good and bad roles in contigs.
print "Analyzing features.\n";
my $count = 0;
my $featureList = $gto->{features};
my $total = scalar(@$featureList);
print "$total features found in GTO.\n";
for my $feature (@$featureList) {
    my $fid = $feature->{id};
    $stats->Add(featureIn => 1);
    # Get this feature's roles.
    my $roles = $fidRoles{$fid};
    if (! $roles) {
        # This feature is not a PEG or does not have an interesting role.
        $stats->Add(featureUnMapped => 1);
    } else {
        # Now get the location.
        my $loc = $feature->{location};
        if (! $loc) {
            # No location, so we don't care.
            $stats->Add(featureNoLocation => 1);
        } else {
            # Are we good, bad, both, or neutral? We go through this rigamarole so that no feature gets counted twice for the same type.
            # So, a feature with 2 bad roles must not be counted as 2 bad features.
            my ($isGood, $isBad) = (0, 0);
            for my $role (@$roles) {
                $stats->Add(roleInFeature => 1);
                if ($ppr{$role}) {
                    $isBad++;
                } else {
                    $isGood++;
                }
            }
            if ($isGood || $isBad) {
                # We have something to record. Get the contig IDs.
                my @contigs = map { $_->[0] } @$loc;
                # Loop through the contigs.
                for my $contig (@contigs) {
                    if (! defined $contigLens{$contig}) {
                        $stats->Add(badContigID => 1);
                    } else {
                        $stats->Add(contigInFeature => 1);
                        if ($isGood) {
                            $contigGood{$contig}++;
                            $stats->Add(goodRoleInContig => 1);
                        }
                        if ($isBad) {
                            push @{$contigBad{$contig}}, $fid;
                            $stats->Add(badRoleInContig => 1);
                        }
                    }
                }
            }
        }
    }
    $count++;
    print "$count of $total features processed.\n" if $count % 1000 == 0;
}
# Now we write the contig report.
open($oh, ">$evalDir/contigs.tbl") || die "Could not open contigs.tbl: $!";
print $oh join("\t", qw(Contig Len Good Bad BadFids)) . "\n";
for my $contig (sort keys %contigLens) {
    my @badFids = @{$contigBad{$contig}};
    print $oh join("\t", $contig, $contigLens{$contig}, $contigGood{$contig}, scalar(@badFids), @badFids) . "\n";
}
close $oh; undef $oh;
print "All done.\n" . $stats->Show();
