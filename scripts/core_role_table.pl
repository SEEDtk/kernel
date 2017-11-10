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
use CopyFromSeed;
use Stats;
use RoleParse;


=head1 Build a Role Table From a FIGdisk

    core_role_table.pl [ options ] coreDir outFile

This script builds a C<roles.in.subsystems>-like file from a FIGdisk's Subsystems directory. All subsystems are included.
(See L<build_role_tables.pl> for a description of the C<roles.in.subsystems> file.) This file can be used to identify a
larger set of subsystem-annotated roles since it will include experimental subsystems. Note, however, that
only roles currently in the L<Shrub> database will be processed, and only exchangeable subsystems will be processed.

The output file will be tab-delimited, with four columns: (0) role ID, (1) role hash code, (2) role name, and (3) "X" if
the role is from an experimental subsystem.

=head2 Parameters

The positional parameters are the name of the FIGdisk directory and the name to give to the output file.

The command-line options are those found in L<Shrub/script_options>.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('coreDir outFile',
            Shrub::script_options(),
        );
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the parameters.
my ($coreDir, $outFile) = @ARGV;
if (! $coreDir) {
    die "No SEED directory specified.";
} elsif (! -d "$coreDir/FIG/Data/Subsystems") {
    die "$coreDir does not appear to be a SEED FIGdisk.";
} elsif (! $outFile) {
    $outFile = "roles.in.all.subsystems";
    print "Default output file name used.\n";
}
# Get the subsystem list.
my $subsysDir = "$coreDir/FIG/Data/Subsystems";
opendir(my $dh, $subsysDir) || die "Could not open subsystem directory: $!";
my @subs = grep { -f "$subsysDir/$_/spreadsheet" } readdir $dh;
closedir $dh;
my $total = scalar @subs;
print "$total subsystem directories found.\n";
# This hash will map each role ID to its checksum and name.
my %roles;
# This hash maps each checksum to a role ID.
my %checksums;
# This will count the subsystems.
my $counter = 0;
# This hash will count the number of non-experimental subsystems containing each role. When we
# unspool the output, roles not in here will get the X at the end.
my %nonExp;
for my $sub (sort @subs) {
    $counter++;
    print "Processing $sub ($counter of $total).\n";
    my $subDir = "$subsysDir/$sub";
    my $public = CopyFromSeed::ReadFlagFile("$subDir/EXCHANGABLE");
    if (! $public) {
        print "Private subsystem skipped.\n";
        $stats->Add(subsystemPrivate => 1);
    } else {
        # This will be set to TRUE if we are experimental.
        my $exp;
        if (open(my $ch, "$subDir/CLASSIFICATION")) {
            my $line = <$ch>;
            if ($line =~ /experimental/i) {
                $exp = 1;
                print "Subsystem $sub is experimental.\n";
                $stats->Add(subsystemExperimental => 1);
            }
        }
        # Now we are ready to read the roles.
        open(my $sh, "$subDir/spreadsheet") || die "Could not open $sub spreadsheet: $!";
        my $line = '';
        while (! eof $sh && $line ne '//') {
            $line = <$sh>;
            chomp $line;
            if ($line =~ /^\S+\t(.+)/) {
                # Here we have a role.
                my $name = $1;
                my $checksum = RoleParse::Checksum($name);
                my $roleID;
                $stats->Add(roleRead => 1);
                if (exists $checksums{$checksum}) {
                    $roleID = $checksums{$checksum};
                    $stats->Add(roleOldFound => 1);
                } else {
                    $stats->Add(roleNewFound => 1);
                    ($roleID) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checksum], 'Role(id)');
                    $checksums{$checksum} = $roleID;
                    if (! $roleID) {
                        $stats->Add(roleNotInShrub => 1);
                    } else {
                        $stats->Add(roleAdded => 1);
                        $roles{$roleID} = [$checksum, $name];
                    }
                }
                # If the role is not experimental, record that fact.
                if ($roleID && ! $exp) {
                    $nonExp{$roleID}++;
                    $stats->Add(roleNotExperimental => 1);
                }
            }
        }
        close $sh;
    }
}
# Now we have processed all the subsystems. Write out the roles.
open(my $oh, ">$outFile") || die "Could not open output file: $!";
for my $roleID (sort keys %roles) {
    my $roleData = $roles{$roleID};
    my $exp = ($nonExp{$roleID} ? '' : 'X');
    print $oh join("\t", $roleID, @$roleData, $exp) . "\n";
    $stats->Add(roleOut => 1);
}
print "All done: " . $stats->Show();