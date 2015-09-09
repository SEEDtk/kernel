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
use Bin;
use Bin::Analyze;
use Shrub;

=head1 Report About Bins

    bins_report.pl [ options ]

Produce a report of the quality bins in a bin list.

=head2 Parameters

The input file is a list of bins,  in JSON format, with the good bins at the top.

The command-line options are those found in L<Shrub::script_options> and L<ScriptUtils/ih_options> plus the following.

=over 4

=item unis

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(),
        Shrub::script_options(),
        ['unis=s',     'universal role file', { default => "$FIG_Config::data/Global/uni_roles.tbl" }],
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
my $binList;
# Read in the universal roles.
open(my $uh, '<', $opt->unis) || die "Could not open universal role file: $!";
my %uniRoles;
while (! eof $uh) {
    my $line = <$uh>;
    chomp $line;
    my ($roleID, undef, $roleName) = split /\t/, $line;
    $uniRoles{$roleID} = $roleName;
}
# Get the analysis object.
my $analyzer = Bin::Analyze->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# This will be a hash mapping each universal role to a hash of the bins it appears in. The bins will be
# identified by an ID number we assign.
my %uniBins;
my $binID = 0;
# This will contain a list of bins for the final statistical report.
my @bins;
# Loop through the bins, producing reports on relatively good ones.
my $done;
while (! eof $ih && ! $done) {
    my $bin = Bin->new_from_json($ih);
    my $quality = $analyzer->AnalyzeBin($bin);
    push @bins, $bin;
    # Compute the bin ID.
    $binID++;
    print "\nBIN $binID (from " . $bin->contig1 . ", " . $bin->contigCount . " contigs, " . $bin->len . " base pairs, quality $quality)\n";
    # List the close reference genomes.
    my @genomes = $bin->refGenomes;
    if (@genomes) {
        my $filter = 'Genome(id) IN (' . join(', ', map { '?' } @genomes) . ')';
        my %gNames =  map { $_->[0] => $_->[1] } $shrub->GetAll('Genome', $filter, \@genomes, 'id name');
        for my $genome (@genomes) {
            print "    $genome: $gNames{$genome}\n";
        }
    }
    # Compute the average coverage.
    my $coverageV = $bin->coverage;
    my $avg = 0;
    for my $covg (@$coverageV) {
        $avg += $covg;
    }
    $avg /= scalar @$coverageV;
    print "*** Mean coverage is $avg.\n";
    # Finally, the universal role list. This hash helps us find the missing ones.
    print "    Universal Roles\n";
    print "    ---------------\n";
    my %unisFound = map { $_ => 0 } keys %uniRoles;
    my $uniFoundCount = 0;
    my $uniMissingCount = 0;
    my $uniDuplCount = 0;
    # Get the universal role hash for the bin.
    my $binUnis = $bin->uniProts;
    for my $uni (sort keys %$binUnis) {
        my $count = $binUnis->{$uni};
        if ($count) {
            print "    $uni\t$uniRoles{$uni}\t$count\n";
            $unisFound{$uni} = 1;
            $uniBins{$uni}{$binID} = $count;
            $uniFoundCount++;
            if ($count > 1) {
                $uniDuplCount++;
            }
        }
    }
    # Now the roles not found.
    if (! $uniFoundCount) {
        print "    NONE FOUND\n";
    } else {
        print "    ---------------\n";
        for my $uni (sort keys %unisFound) {
            if (! $unisFound{$uni}) {
                print "    $uni\t$uniRoles{$uni}\tmissing\n";
                $uniMissingCount++;
            }
        }
        print "    ---------------\n";
        print "    $uniFoundCount present, $uniMissingCount missing, $uniDuplCount duplicated.\n";
    }
}
# Now output the universal role matrix.
print "\nUNIVERSAL ROLE MATRIX\n";
print join("\t", 'Role', map { "bin$_" } (1 .. $binID)) . "\n";
for my $uni (sort keys %uniBins) {
    print join("\t", $uni, map { $uniBins{$uni}{$_} // ' ' } (1 .. $binID)) . "\n";
}
print "\n\n";
# Finally, the bin statistics.
my $stats = $analyzer->Stats(\@bins);
print "FINAL REPORT\n\n" . $stats->Show();