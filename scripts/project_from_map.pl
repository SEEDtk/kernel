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
use MapToRef;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use JSON::XS;

=head1 project reference genome to a close strain

    project_from_map SkeletalGenomeDir [ options ]

project a reference genome to call features


=head2 Parameters

## describe positional parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item -g SkeletalGenomeDir

a path to a skeletal SEED genome directory that must include

    map.json from a previous build_ref_map run

=back

=cut

my $k = 30;    # kmers for generating map.  Chosen conservatively.

# Get the command-line parameters.
my $opt = ScriptUtils::Opts(
	'',
	[ 'genome|g=s',    'Path to a skeletal genome directory' ]
);
my $genomeD = $opt->genome;

my $json = JSON::XS->new;
my $in_fh;

my $refD;

if (open($in_fh,  "<", "$genomeD/ref_dir")) {
    local $/;
    undef $/;
    $refD = <$in_fh>;
    chop $refD;
    close ($in_fh);
} else {
    die "Cannot open $genomeD/ref_dir";
}

my $genetic_code = &get_genetic_code($genomeD);

my $map;
if (open($in_fh, "<", "$genomeD/map.json")) {
    local $/;
    undef $/;
    my $map_txt = <$in_fh>;
    $map = $json->decode($map_txt);
    close ($in_fh);
} else {
    die "Cannot open $genomeD/map.json";
}

if (0) {
        foreach my $m (@$map) {
            my @tmp = (@{$m->[0]}, @{$m->[1]});
            print join("\t", @tmp), "\n";
        }
}
my @g_tuples;
if (open($in_fh, "<", "$genomeD/gtuples")) {
    local $/; 
    undef $/; 
    my $g_txt = <$in_fh>;
    @g_tuples = $json->decode($g_txt);
    close($in_fh);
} else {
        die "Cannot open $genomeD/gtuples";
}

&MapToRef::build_features( $map, $refD, $genomeD, \@g_tuples, $genetic_code );


sub get_genetic_code {
	my ($dir) = @_;

	if ( !-s "$dir/GENETIC_CODE" ) { return 11 }
	my @tmp = `cat $dir/GENETIC_CODE`;
	chomp $tmp[0];
	return $tmp[0];
}
