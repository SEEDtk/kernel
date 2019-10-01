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

=head1 Collate Features by Distance

    p3x-close-check.pl [options] pegCol roleCol

This script takes an input file containing feature IDs and role IDs.  It compares the features against the database to determine
how close they are.  Features within a certain distance of each other are written to the output.  Input lines with a blank role ID
will be ignored.

=head2 Parameters

The positional parameters are the names/indices of the feature ID column and the role ID column.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  It should be tab-delimited.

Additional command-line options are those given in L<Shrub/script_options> to select the L<Shrub> database
plus the following.

=over 4

=item verbose

Display progress messages on STDERR.

=item gap

The maximum gap distance for two features to be considered neighbors.  The default is 10000.

=cut

use strict;
use P3Utils;
use Shrub;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('pegCol roleCol', Shrub::script_options(), P3Utils::ih_options(),
        ['verbose|debug|v', 'display progress messages on STDERR'],
        ['gap|g=i', 'maximum gap distance between neighbors', { default => 10000 }]
);

# Get access to the Shrub.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = P3Utils::ih($opt);
my ($pegCol, $roleCol) = @ARGV;
if (! defined $pegCol || ! defined $roleCol) {
    die "ID and role columns are required parameters.";
}
# Get the options.
my $debug = $opt->verbose;
my $gap = $opt->gap;
# Read the incoming headers.
my ($headers, $cols) = P3Utils::find_headers($ih, input => $pegCol, $roleCol);
# Write the output headers.
P3Utils::print_cols(['id1', 'id2', 'distance', 'role1', 'role2']);
# This hash will track the role descriptions.
my %roles;
# This hash will map each feature ID to a role.
my %pegRoles;
# This hash will map each feature ID to a location.
my %pegLocs;
# Loop through the input.
print STDERR "Reading input features.\n" if $debug;
while (! eof $ih) {
    my ($fid, $roleID) = P3Utils::get_cols($ih, $cols);
    # Only process features with a role.
    if ($roleID) {
        print STDERR "$fid found with role $roleID.\n" if $debug;
        # Get the role description.
        my $role = $roles{$roleID};
        if (! $role) {
            $role = $shrub->role_id_to_desc($roleID);
            $roles{$roleID} = $role;
        }
        $pegRoles{$fid} = $role;
        # Get the location.
        $pegLocs{$fid} = $shrub->loc_of($fid);
    }
}
# Now find the neighboring features.
my @fids = sort keys %pegRoles;
print STDERR scalar(@fids) . " features found with " . scalar(keys %roles) . " roles.\n" if $debug;
while (@fids) {
    my $fid1 = shift @fids;
    my $loc1 = $pegLocs{$fid1};
    my $role1 = $pegRoles{$fid1};
    for my $fid2 (@fids) {
        my $distance = $loc1->Distance($pegLocs{$fid2});
        if (defined $distance && $distance < $gap) {
            P3Utils::print_cols([$fid1, $fid2, $distance, $role1, $pegRoles{$fid2}]);
        }
    }
}
print STDERR "All done.\n" if $debug;
