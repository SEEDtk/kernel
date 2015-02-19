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

=head1 Build FIGfams from privileged Subsystems

    build_ff_from_reliable_ss [options] > fam.peg.func

This script builds FIGfams from privileged Subsystems.  It
assumes that all PEGs in a column from a subsystem that have
identical privileged Functions are, in fact, homologous.  That is, it
assumes that curators use distinct Roles for convergent versions
of an "abstract role".

Each family will contain isofunctional PEGs (that is, they have the same
privileged Function, and they have no comment).  Further each PEG in a
family occurs in at least one row in a populated subsystem that has an
"active variant code (i.e., not -1 or '*-1').

To make a normalize set of families use

    build_ff_from_reliable_ss > tmp
    cut -f1,2 tmp > families.2c
    cut -f1,3 tmp | sort -u > family.functions
    rm tmp

The families.2c and family.functions files are the common form of
a FIGfam release.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options>.

=cut

use strict;
use Data::Dumper;
use Shrub;
use ScriptUtils;

# Get the command-line options.
my $opt = ScriptUtils::Opts('', Shrub::script_options());
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get a list of all the privileged subsystems.
my @privileged_subsystems = map { $_->[0] }
                            $shrub->GetAll("Subsystem","Subsystem(security) = ?",[Shrub::PRIV],"Subsystem(id)");
my %poss_fams;
foreach my $ss (@privileged_subsystems)
{
    # Get all the features in the subsystem. Note that vacant subsystems are not included in the database,
    # so we don't need to filter on variant codes any more.
    my $fidList = $shrub->Subsystem2Feature(id => $ss);
    # Get the related functions.
    my $funHash = $shrub->Feature2Function(Shrub::PRIV, $fidList);
    foreach my $peg ($fidList)
    {
        my ($funID, $function, $comment) = @{$funHash->{$peg}};
        # Discard features with "trunc" or "frame" in the comment.
        unless ($comment =~ /trunc|frame/) {
            push(@{$poss_fams{$function}},$peg);
        }
    }
}
my $nxt_fam = "FIG00000001";
# Get the function names sorted by the number of pegs in each function's group,
# biggest group first.
my @functions = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [$_,scalar @{$poss_fams{$_}}] } keys(%poss_fams);
foreach my $f (@functions)
{
    # Get the list of pegs in this family.
    my $pegsL = $poss_fams{$f};
    if (scalar(@$pegsL) > 1)
    {
        foreach my $peg (sort { &SeedUtils::by_fig_id($a,$b) } (@$pegsL))
        {
            print join("\t",($nxt_fam,$peg,$f)),"\n";
        }
        $nxt_fam++;
    }
}

