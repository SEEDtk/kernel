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

=head1 Get estimates of number of bins for each universal roles

    tabulate_bin_estimates [ options ] < contigs.data

For each universal role, you have an estimate of the number of bins
given by the number of contigs with that role.  The assumption is that
all bins for contigs from genomes with a coverage greater than 10 will include 
one hit for each universal role (which is not completely true).

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options());
my $ih = ScriptUtils::IH($opt->input);

use Data::Dumper;
use Bin;

my $contigs = &Bin::ReadContigs($ih);

my %counts;
foreach my $contig (@$contigs)
{
    my $uni_protH = $contig->uniProts;
    foreach my $role (keys(%$uni_protH))
    {
	$counts{$role}++;
    }
}
use Shrub;
my $shrub = Shrub->new();
foreach my $role (sort { $counts{$b} <=> $counts{$a} } keys(%counts))
{
    my $desc = $shrub->role_id_to_desc($role);
    print join("\t",($counts{$role},$role,$desc)),"\n";
}
