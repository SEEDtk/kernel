#!/usr/bin/env perl

#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
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

The command-line options are those found in L<Shrub/new_for_script>.

=cut

use strict;
use Data::Dumper;
use Shrub;

# Connect to the database.
my ($shrub, $opt) = Shrub->new_for_script('%c %o', { });
my @privileged_subsystems = map { $_->[0] }
                            $shrub->GetAll("Subsystem","Subsystem(security) = ?",[2],"Subsystem(id)");
foreach my $ss (@privileged_subsystems)
{
    my %active_genomes = map { ($_->[0] => 1) }
                         $shrub->GetAll("Subsystem2Genome",
          "(Subsystem2Genome(from-link) = ? AND Subsystem2Genome(variant) != ?) and (Subsystem2Genome(variant) != ?)",
          [$ss, '-1', '*-1'],
          "Subsystem2Genome(to-link)");
    my %poss_fams;
    ## NOTE: Feature2Function does not always work, because most functions connect to proteins. I created a Shrub method
    ##       to circumvent this.
    # Get all the features in the subsystem.
    my @fids = map { $_->[0] } $shrub->GetAll("Subsystem2Feature", 'Subsystem2Feature(from-link) = ?', [$ss], 'Subsystem2Feature(to-link)');
    # Get the related functions.
    my $funHash = $shrub->Feature2Function(2, \@fids);
    foreach my $peg (@fids)
    {
        my $function = $funHash->{$peg}[1];
        # Discard features with "trunc" or "frame" in the comment.
        unless ($function =~ /trunc|frame/) {
            my $g = &SeedUtils::genome_of($peg);
            if ($active_genomes{$g})
            {
                push(@{$poss_fams{$function}},$peg);
            }
        }
    }
    my $nxt_fam = "FIG00000001";
    my @functions = map { $_->[0] }
                    sort { $b->[1] <=> $a->[1] }
                    map { [$_,scalar @{$poss_fams{$_}}] } keys(%poss_fams);
    foreach my $f (@functions)
    {
        my $pegsH = $poss_fams{$f};
        if (keys(%$pegsH) > 1)
        {
            foreach my $peg (sort { &SeedUtils::by_fig_id($a,$b) } keys(%$pegsH))
            {
                print join("\t",($nxt_fam,$peg,$f)),"\n";
            }
            $nxt_fam++;
        }
    }
}
