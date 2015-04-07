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
use Data::Dumper;
use JSON::XS;
use Shrub;
use ScriptUtils;
use ClosestReps;

my $k = 10;

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '', 
                     Shrub::script_options(), ScriptUtils::ih_options(),
                        []
    );
my $ih = ScriptUtils::IH( $opt->input );
my $shrub = Shrub->new_for_script($opt);

my %g_to_gs = map { ($_->[0] => $_->[1]) }
              map { chomp; [split(/\t/,$_)] }
              <$ih>;

my @roles = <DATA>;
chomp @roles;
my %roles = map { ($_ => 1) } @roles;
my $kmer_hash = &Closest_reps::load_kmers( \%roles, \%g_to_gs, $k, $shrub );
my $close_ref_data = { kmers => $kmer_hash, genomes => \%g_to_gs };
&SeedUtils::write_encoded_object($close_ref_data,\*STDOUT);


__DATA__
Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
Phenylalanyl-tRNA synthetase beta chain (EC 6.1.1.20)
Preprotein translocase secY subunit (TC 3.A.5.1.1)
GTP-binding and nucleic acid-binding protein YchF
Translation initiation factor 2
Signal recognition particle, subunit Ffh SRP54 (TC 3.A.5.1.1)
Histidyl-tRNA synthetase (EC 6.1.1.21)
Methionyl-tRNA synthetase (EC 6.1.1.10)
Isoleucyl-tRNA synthetase (EC 6.1.1.5)
Valyl-tRNA synthetase (EC 6.1.1.9)
Seryl-tRNA synthetase (EC 6.1.1.11)
Alanyl-tRNA synthetase (EC 6.1.1.7)
