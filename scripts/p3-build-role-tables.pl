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
use RoleParse;
use CopyFromSeed;
use Stats;
use SeedUtils;

=head1 Build Role Tables From the Patric Snapshot

    p3-build-role-tables.pl [ options ] annotations subsystems outDir

This script builds the same role tables output by L<build_role_tables.pl>, but builds them from the PATRIC
kmer snapshot rather than the Shrub. The goal is to create something closer to the data used by RASTtk.

PATRIC does not have a real snapshot of the coreSEED. Instead, there is a compressed copy of the subsystem
directory and a directory containing assignments for each genome. The subsystem directory must be
decompressed before this script is used.

The assignment directory contains a file for each genome, with the same name as the genome ID. The file is
tab-delimited. Each row contains (0) a feature ID and (1) a functional assignment.

The Shrub database is used to convert checksums to role IDs.

Note that the C<roles.to.use> file will be restricted to non-hypothetical roles that occur between 1 to 5 times in over 100
genomes.



=head2 Parameters

The positional parameters are the name of the assignment directory, the name of the subsystem directory, and
the name of the output directory.

Annotations are in C</vol/core-seed/kmers/core.>I<XXXX-XXXX>C</Annotations/0>. Subsystems are in
C</vol/patric3/fams/>I<XXXX-XXXX>C</subsystem-import/subsystems.tgz>. In both cases I<XXXX-XXXX> is
the release date. The subsystems need to be decompressed before processing.

The command-line options are those found in L<Shrub/script_options>.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('annotations subsystems outDir',
        Shrub::script_options(),
        );
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get a list of acceptable genomes.
my %genomes = map { $_ => 1 } $shrub->GetFlat('Genome', 'Genome(well-behaved) = ?', [1], 'id');
# Get the directories.
my ($annotations, $subsystems, $outDir) = @ARGV;
if (! $annotations) {
    die "No annotation directory specified.";
} elsif (! -d $annotations) {
    die "Annotation directory $annotations not found.";
} elsif (! $subsystems) {
    die "No subsystem directory specified.";
} elsif (! -d $subsystems) {
    die "Subsystem directory $subsystems not found.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Output directory $outDir not found.";
}
# Get the list of subsystem directories.
print "Processing subsystems in 4subsystems.\n";
opendir(my $sdh, $subsystems) || die "Could not open $subsystems: $!";
my @subsystems = grep { -s "$subsystems/$_/spreadsheet " } readdir $sdh;
closedir $sdh;
print scalar(@subsystems) . " subsystems found.\n";
# This hash will contain the roles found in subsystems. For each checksum, it will contain the role ID and name.
my %roleHash;
# Loop through the subsystems.
for my $subsystem (@subsystems) {
    my $dirName = "$subsystems/$subsystem";
    if (! CopyFromSeed::ReadFlagFile("$dirName/EXCHANGABLE")) {
        print "Subsystem $subsystem is private in SEED.\n";
        $stats->Add(subsystemPrivate => 1);
    } else {
        # Check the classification for an experimental subsystem.
        my $ok = 1;
        if (open(my $ch, "$dirName/CLASSIFICATION")) {
            my $line = <$ch>;
            if ($line =~ /experimental/i) {
                $ok = 0;
                print "Subsystem $subsystem is experimental.\n";
                $stats->Add(subsystemExperimental => 1);
            }
        }
        if ($ok) {
            # This is a real subsystem. Process it.
            print "Reading roles from $subsystem.\n";
            open(my $sh, "<$dirName/spreadsheet") || die "Could not open $subsystem spreadsheet: $!";
            my $line = <$sh>;
            while ($line && $line =~ /^\S+\t(.+)/) {
                my $function = $1;
                $stats->Add(roleLineIn => 1);
                my @roleNames = SeedUtils::roles_of_function($function);
                for my $roleName (@roleNames) {
                    my $checksum = RoleParse::Checksum($roleName);
                    if ($roleHash{$checksum}) {
                        # Here we have already stored this role.
                        $stats->Add(roleDuplicate => 1);
                    } else {
                        my ($roleID) = $shrub->GetFlat('Role', 'Role(checksum) = ? AND Role(hypo) = ?', [$checksum, 0], 'id');
                        if(! $roleID) {
                            # Here we have a role that is no longer supported, which is very bad.
                            print "WARNING: No role ID found for $roleName.\n";
                            $stats->Add(roleNotFound => 1);
                        } else {
                            $roleHash{$checksum} = [$roleID, $roleName];
                            $stats->Add(roleStored => 1);
                        }
                    }
                }
                $line = <$sh>;
            }
            $stats->Add(subsystemProcessed => 1);
        }
    }
}
my $roleCount = scalar keys %roleHash;
if ($roleCount == 0) {
    die "No roles found.\n";
} else {
    print "$roleCount roles found.\n";
}
# All the subsystem roles are done. Build roles.in.subsystems.
open(my $oh, ">$outDir/roles.in.subsystems") || die "Could not open roles.in.subsystems: $!";
for my $checksum (sort keys %roleHash) {
    my ($roleID, $roleName) = @{$roleHash{$checksum}};
    print $oh join("\t", $roleID, $checksum, $roleName) . "\n";
}
close $oh; undef $oh;
# Now we build the raw table. Get the genomes.
opendir(my $adh, $annotations) || die "Could not open annotations directory: $!";
my @genomes = sort grep { $genomes{$_} && -s "$annotations/$_" } readdir $adh;
closedir $adh;
# Open the output file.
open($oh, ">$outDir/raw.table") || die "Could not open raw.table: $!";
# The role counts will be stored in here.
my %roleCount;
# Loop through the genomes.
for my $genome (@genomes) {
    print "Processing $genome.\n";
    open(my $ih, "$annotations/$genome") || die "Could not open $genome annotations: $!";
    # This will track the genome counts for each role.
    my %gRoleCount;
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(genomeLineIn => 1);
        if ($line =~ /^(\S+)\t(.+)/) {
            my ($peg, $function) = ($1, $2);
            my @roleNames = SeedUtils::roles_of_function($function);
            for my $roleName (@roleNames) {
                my $checksum = RoleParse::Checksum($roleName);
                if (! exists $roleHash{$checksum}) {
                    $stats->Add(annotationRoleSkipped => 1);
                } else {
                    my $roleID = $roleHash{$checksum}[0];
                    $gRoleCount{$roleID}++;
                    print $oh join("\t", $roleName, $roleID, $peg) . "\n";
                    $stats->Add(annotationRoleOut => 1);
                }
            }
        }
    }
    # Loop through the roles in this genome, counting the ones that occur less than 5 times.
    for my $roleID (keys %gRoleCount) {
        my $rCount = $gRoleCount{$roleID};
        if ($rCount >= 1 && $rCount <= 5) {
            $roleCount{$roleID}++;
        }
    }
}
close $oh; undef $oh;
# Now output the roles that occur 100 times or more.
print "Writing roles.to.use.\n";
open($oh, ">$outDir/roles.to.use") || die "Could not open roles.to.use: $!";
for my $roleID (sort keys %roleCount) {
    if ($roleCount{$roleID} >= 100) {
        print $oh "$roleID\n";
        $stats->Add(roleAcceptedForUse => 1);
    }
}
print "All done.\n" . $stats->Show();