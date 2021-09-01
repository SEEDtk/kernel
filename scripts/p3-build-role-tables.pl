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
use P3DataAPI;
use P3Utils;

=head1 Build Role Tables From the Patric Snapshot

    p3-build-role-tables.pl [ options ] outDir annotations subsystems

This script builds the same role tables output by L<build_role_tables.pl>, but builds them from the BV-BRC
kmer snapshot rather than the Shrub. The goal is to create something closer to the data used by RASTtk.

BV-BRC does not have a real snapshot of the coreSEED. Instead, there is a compressed copy of the subsystem
directory and a directory containing assignments for each genome. The subsystem directory must be
decompressed before this script is used.

The assignment subdirectories contains a file for each genome, with the same name as the genome ID. The file is
tab-delimited. Each row contains (0) a feature ID and (1) a functional assignment.

The Shrub database is used to convert checksums to role IDs.

Note that the C<roles.to.use> file will be restricted to non-hypothetical roles that occur between 1 to 5 times in over 100
genomes.

=head2 Basic Procedure for Creating the Predictors.

First, you must identify a directory for the predictors and a directory for the role files. In the example below,
we will use C<FunctionPredictors.X> for the predictors and C<Staging> for the role files.

The first command creates the initial role files.

    p3-build-role-tables Eval.XX Annotations Subsystems

The following three commands are run repeatedly in a loop. When the output of L<build_roles_to_use.pl> indicates that
all of the roles are good, the loop stops.

    build_matrix --clear Eval.XX Eval.XX/FunctionPredictors
    build_LDA_predictors Eval.XX/FunctionPredictors
    build_roles_to_use Eval.XX/FunctionPredictors Eval.XX

Once the role set has converged, to install the new predictors you must (1) copy C<roles.in.subsystems> from C<Staging>
to the SEEDtk P3Data directory. (2) Copy C<FunctionPredictors.X/roles.to.use> to the SEEDtk P3Data directory. (3)
Symbolically link C<FunctionPredictors> in the P3Data directory to C<FunctionPredictors.X>.

=cut

=head2 Parameters

The positional parameters are the name of the output directory, the name of the assignment directory, and
the name of the subsystem directory.

Annotations are in C</vol/core-seed/kmers/core.>I<XXXX-XXXX>C</Annotations/0>. Subsystems are in
C</vol/patric3/fams/>I<XXXX-XXXX>C</subsystem-import/subsystems.tgz>. In both cases I<XXXX-XXXX> is
the release date. The subsystems need to be decompressed before processing. If the subsystem directory
is omitted, then it is assumed C<roles.in.subsystems> already exists in the output directory.  If both
directories are omitted, then it is assumed all genomes are being loaded from BV-BRC.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item patric

The name of a file containing BV-BRC genome IDs in the first column.  These will be processed in addition to
the coreSEED genomes.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir annotations subsystems',
        Shrub::script_options(),
        ['patric=s', 'name of a file containing BV-BRC genome IDs'],
        );
my $stats = Stats->new();
# Get the directories.
my $core = 1;
my ($outDir, $annotations, $subsystems) = @ARGV;
if (! $annotations) {
    $core = 0;
} elsif (! -d $annotations) {
    die "Annotation directory $annotations not found.";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Output directory $outDir not found.";
}
if (! $subsystems) {
    if (! -s "$outDir/roles.in.subsystems") {
        die "No subsystem directory specified and no roles.in.subsystems present.";
    }
} elsif (! -d $subsystems) {
    die "Subsystem directory $subsystems not found.";
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get a list of acceptable genomes.
my %genomes;
if ($core) {
    %genomes = map { $_ => 1 } $shrub->GetFlat('Genome', 'Genome(well-behaved) = ?', [1], 'id');
}
# This hash will contain the roles found in subsystems. For each checksum, it will contain the role ID and name.
my %roleHash;
# Get the list of subsystem directories.
if ($subsystems) {
    print "Processing subsystems in $subsystems.\n";
    opendir(my $sdh, $subsystems) || die "Could not open $subsystems: $!";
    my @subsystems = grep { -s "$subsystems/$_/spreadsheet" } readdir $sdh;
    closedir $sdh;
    print scalar(@subsystems) . " subsystems found.\n";
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
} else {
    # Here we read the roles from an existing roles.in.subsystems.
    open(my $ih, "<$outDir/roles.in.subsystems") || die "Could not open roles.in.subsystems: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\S+)\t(\S+)\t(.+)/) {
            $roleHash{$2} = [$1, $3];
        }
    }
    print scalar(keys %roleHash) . " roles read from file.\n";
}
# Now we build the raw table. Get the genomes.
my %genomeFiles;
if ($core) {
    opendir(my $dh, $annotations) || die "Could not open annotations directory: $!";
    while (my $file = readdir $dh) {
        if ($genomes{$file}) {
            my $path = "$annotations/$file";
            if (-s $path) {
                $genomeFiles{$file} = $path;
            }
        }
    }
}
# Open the output file.
open(my $oh, ">$outDir/raw.table") || die "Could not open raw.table: $!";
# The role counts will be stored in here.
my %roleCount;
# The bad roles will be listed in here.
my %badRoles;
# Loop through the genomes in files.
for my $genome (sort keys %genomeFiles) {
    print "Processing $genome from CoreSEED.\n";
    open(my $ih, '<', $genomeFiles{$genome}) || die "Could not open $genome annotations: $!";
    # This will track the genome counts for each role.
    my %gRoleCount;
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(featureLineIn => 1);
        if ($line =~ /^(\S+)\t(.+)/) {
            my ($peg, $function) = ($1, $2);
            processRoles($peg, $function, \%gRoleCount);
        }
    }
    close $ih; undef $ih;
    countRoles(\%gRoleCount);
}
# Loop through the genomes in BV-BRC.
if ($opt->patric) {
    my $p3 = P3DataAPI->new();
    open(my $ih, '<', $opt->patric) || die "Could not open BV-BRC genome file: $!";
    my @genomes;
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\d+\.\d+)/) {
            my $genome = $1;
            my %gRoleCount;
            print "Processing $genome from BV-BRC.\n";
            my $features = P3Utils::get_data($p3, feature => [['eq', 'genome_id', $genome]], ['patric_id', 'product']);
            for my $feature (@$features) {
                my ($peg, $function) = @$feature;
                processRoles($peg, $function, \%gRoleCount);
            }
            countRoles(\%gRoleCount);
        }
    }
    close $ih; undef $ih;
}
close $oh; undef $oh;
# Now output the roles that occur 100 times or more.
print "Writing roles.to.use.\n";
open($oh, ">$outDir/roles.to.use") || die "Could not open roles.to.use: $!";
for my $roleID (sort keys %roleCount) {
    if ($roleCount{$roleID} >= 100 && ! $badRoles{$roleID}) {
        print $oh "$roleID\n";
        $stats->Add(roleAcceptedForUse => 1);
    }
}
print "All done.\n" . $stats->Show();

# Count the roles in a peg and write them to the output.
sub processRoles {
    my ($peg, $function, $gRoleCount) = @_;
    my @roleNames = SeedUtils::roles_of_function($function);
    for my $roleName (@roleNames) {
        my $checksum = RoleParse::Checksum($roleName);
        if (! exists $roleHash{$checksum}) {
            $stats->Add(annotationRoleSkipped => 1);
        } else {
            my $roleID = $roleHash{$checksum}[0];
            $gRoleCount->{$roleID}++;
            print $oh join("\t", $roleName, $roleID, $peg) . "\n";
            $stats->Add(annotationRoleOut => 1);
        }
    }
}

# Fold a genome's roles into the counter.
sub countRoles {
    my ($gRoleCount) = @_;
    # Loop through the roles in this genome, keeping the ones that occur 5 times or less.
    for my $roleID (keys %$gRoleCount) {
        my $rCount = $gRoleCount->{$roleID};
        if ($rCount > 5) {
            $badRoles{$roleID} = 1;
        } elsif ($rCount >= 1) {
            $roleCount{$roleID}++;
        }
    }

}
