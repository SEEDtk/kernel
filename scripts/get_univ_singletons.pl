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
use JSON::XS;
use warnings;
use Shrub;
use ScriptUtils;
use Data::Dumper;

=head1 Get Roles that are universal and occur as singletons in prokaryotic genomes

    get_univ_singletons [ options ] > RolesAndCounts

=head2 Parameters

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '', Shrub::script_options(),
    [] );

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my @tuples = $shrub->GetAll("Function Function2Feature",
                            "Function2Feature(security) = ?",[1],
                            "Function(description) Function2Feature(to-link) Function2Feature(from-link)");

my %funcH;
my %funcD;
foreach my $tuple (@tuples)
{
    my($funcD,$fid, $func) = @$tuple;
    if (! &SeedUtils::hypo($funcD))
    {
        if ($fid =~ /^fig\|(\d+\.\d+)\.peg\./)
        {
            $funcH{$func}->{$1}++;
            $funcD{$func} = $funcD;
        }
    }
}
my (%counts, %multis);
foreach my $func (keys(%funcH))
{
    my $genomesH = $funcH{$func};
    my @genomes = keys(%$genomesH);
    foreach my $g (keys(%$genomesH))
    {
        my $occur = $genomesH->{$g};
        if ($occur == 1) {
            $counts{$func}++
        } elsif ($occur > 1) {
            $multis{$func}++;
        }
    }
}

my @best = sort { $counts{$b} <=> $counts{$a} } grep { ($multis{$_} // 0)  < 10 } keys(%counts);
foreach $_ (@best)
{
    print join("\t",($_, $counts{$_}, $funcD{$_})),"\n";
}
