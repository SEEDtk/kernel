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
use Data::Dumper;
use Shrub;
use ScriptUtils;
use RolesInSubsystems;

=head1 Create Coupling Report

    non_ss_coupled_to_subsystem.pl [ options ]

This script creates the final coupling report. It searches for roles not in subsystems that are coupled to subsystem roles.
The output file is tab-delimited, and divided into sections. Each section contains a header line specifying the non-subsystem
role and one or more data lines specifying the clustered roles.

The header line fields are (0) the maximum coupling count for
the role, (1) the role description, and (2) the ID of a feature that exemplifies the role. This exemplar is chosen from features
that are nearby as many coupled roles as possible.

The data line fields are (0) a blank, indicating that this is a data line, (1) the name of the coupled role, and (2) an C<X> if the role
is in an experimental subsystem.

A line containing a separator (C<//>) divides sections with exemplar features from those without.

=head2 Parameters

There are no positional parameters.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item gap

Region size in each direction for finding clustered features. The default is C<5000>.

=item roles

The input file specifying which roles are in subsystems. The default is C<roles.in.all.subsystems>. This file is tab-delimited, and
each record contains (0) a role ID, (1) a role checksum, (2) the role name, and (3) an C<X> if the role is in an experimental subsystem.
The default is C<roles.in.all.subsystems>.

=item roleCouples

A file containing the role couples. This file is tab-delimited with headers, and each record contains (0) the first role name, (1) the
second role name, and (2) the number of coupling instances. The default is C<role.couples.tbl>.

=item min

The minimum number of coupling instances required for a pair to be significant. The default is C<100>.

=back

=head2 Procedure

To generate the input files for this program, create a special directory and use the following commands.

    all_entities Genome | get_feature_coupling_table –-headers >coupling.tbl 2>coupling.log
    p3-find-couples --location=location --sequence=contig role <coupling.tbl >role0.couples.tbl
    expand_coupled_roles <role0.couples.tbl >role.couples.tbl
    core_role_table /vol/core-seed/FIGdisk roles.in.all.subsystems >roles.log
    non_ss_coupled_to_subsystem >coupling.report.txt

=cut


$| = 1;
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        ['gap=i', 'region size in each direction', { default => 5000 }],
        ['roles=s', 'roles.in.subsystems file', { default => 'roles.in.all.subsystems'}],
        ['roleCouples=s', 'role couples file', { default => 'role.couples.tbl' }],
        ['min=i', 'minimum count', { default => 100 }]);
my $roles_in_ssF  = $opt->roles;
my $role_couplesF = $opt->rolecouples;
my $min_count     = $opt->min;
my $shrub         = Shrub->new_for_script($opt);
my $gap           = $opt->gap;
print STDERR "Creating role map.\n";
my $roleMap = RolesInSubsystems->new($shrub, $roles_in_ssF);

my %poss;
print STDERR "Reading couples file.\n";
my $progress = 0;
open(my $rh, "<$role_couplesF") || die "could not open $role_couplesF";
# Skip header line.
my $line = <$rh>;
# Loop through the couples.
while (defined($line = <$rh>)) {
    chomp $line;
    my($r1,$r2,$count) = split(/\t/,$line);
    next if ($count < $min_count);
    next if (&ignore($r1) || &ignore($r2));
    my $id1 = $roleMap->RoleID($r1);
    my $id2 = $roleMap->RoleID($r2);
    $progress++;
    print STDERR "$progress couples processed.\n" if ($progress % 1000) == 0;
    # Insure both roles are valid (basically, they are not hypothetical or ill-formed).
    if ($id1 && $id2) {
        my $sub1 = $roleMap->in_sub($id1);
        my $sub2 = $roleMap->in_sub($id2);
        if ($sub1 && ! $sub2) {
            $poss{$id2}{$id1} = $count;
        } elsif (! $sub1 && $sub2) {
            $poss{$id1}{$id2} = $count;
        }
    }
}
close $rh; undef $rh;
# %poss{r1}{r2} now contains the number of couples of non-subsystem role r1 with subsystem role r2.
my @tocheck;
print STDERR "Creating check queue.\n";
foreach my $r (sort keys(%poss))
{
    push(@tocheck,[$r,$poss{$r},&best_count($poss{$r})]);
}

my @sorted = sort { $b->[2] <=> $a->[2] } @tocheck;
my @report2;
my $reportL;
print STDERR scalar(@sorted) . " couple groups to check.\n";
$progress = 0;
foreach my $tuple (@sorted)
{
    my($r1,$coupled,$count) = @$tuple;
    $progress++;
    print "$progress couple groups checked.\n" if ($progress % 100) == 0;
    my $peg = &exemplar($r1,$coupled,$roleMap);
    if ($peg) {
        $reportL = [];
    } else {
        $reportL = \@report2;
    }
    push @$reportL, join("\t", $count, $roleMap->name($r1), $peg) . "\n";
    my @couples = sort { $coupled->{$b} <=> $coupled->{$a} } keys %$coupled;
    for my $couple (@couples) {
        push @$reportL, join("\t", $coupled->{$couple}, $roleMap->name($couple), $roleMap->is_exp($couple)) . "\n";
    }
    if ($peg) {
        print @$reportL;
    }
}
print  "//", @report2;


sub exemplar {
    my($r1,$coupled,$roleMap) = @_;
    # The use of "2" restricts us to CoreSEED functions.
    my %pegs = map { $_ => 0 } $shrub->GetFlat('Function2Feature', 'Function2Feature(from-link) = ? AND Function2Feature(security) = ?',
            [$r1, 2], 'Function2Feature(to-link)');
    # Get the IDs of all the roles in which we are interested.
    my %roles = map { $roleMap->RoleID($_) => 1 } keys %$coupled;
    for my $peg (keys %pegs) {
        # Get the location of this peg.
        my $loc = $shrub->loc_of($peg);
        # Expand it by the gap distance. We might fall off the end of the contig, but we don't care, since there are no pegs there.
        $loc->Widen($gap);
        # Get all the genes in the region. The function ID is in list position 2 of each tuple.
        my $geneList = $shrub->genes_in_region($loc, 2);
        # Loop through the genes counting roles.
        for my $geneThing (@$geneList) {
            for my $role (split /[\@\;\/]/, $geneThing->[2]) {
                if ($roles{$role}) {
                    $pegs{$peg}++;
                }
            }

        }
    }
    # Find the best peg.
    my ($retVal, $count) = ('', 0);
    for my $peg (keys %pegs) {
        if ($pegs{$peg} > $count) {
            $retVal = $peg;
            $count = $pegs{$peg};
        }
    }
    return $retVal;
}

sub best_count {
    my($coupledH) = @_;

    my $best_sofar = 0;
    foreach my $close (keys(%$coupledH))
    {
        if ((! $best_sofar) || ($coupledH->{$close} > $best_sofar))
        {
            $best_sofar = $coupledH->{$close};
        }
    }
    return $best_sofar;

}

sub ignore {
    my($r) = @_;

    if ($r =~ /transport|permease|family/i)
    {
        return 1;
    }
    if ($r =~ /mobile element/i)
    {
        return 1;
    }
    return 0;
}
