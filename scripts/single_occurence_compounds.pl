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
use ScriptUtils;
use Shrub;

=head1 Given a Metobolic-Network (a set of reactions), Compute Compounds that Occur in just one Reaction

Given a putative metabolic network, look for compounds that occur in just
one Reaction.  These represent spots where something is probably amiss.

=head2 Parameters

The standard input must be tab-delimited and contain Reaction IDs in one column.  Only Reactions
that contain at least one single-ocurrence compound will be retained, and a new column
containing a singly-occurring compound will be added as the last field in each retained row.

The input file is either the standard input or specified via L<ScriptUtils::ih_options>.

The command-line options are those in L<Shrub/script_options> plus the following.

=over 4

=item col

Index (1-based) of the column containing the reaction IDs. The default is C<0>, indicating the
last column.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('', ["col|c=i", 'input column index', { default => 0 }],
        Shrub::script_options(), ScriptUtils::ih_options());
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
my %reactions_containing;
# Get the reaction IDs.
my @couplets = ScriptUtils::get_couplets($ih, $opt->col, 0);
for my $reaction (map { $_->[0] } @couplets) {
    my @compounds = $shrub->GetFlat('Reaction2Compound', 'Reaction2Compound(from-link) = ?',
            [$reaction], 'to-link');
   for my $compound (@compounds) {
       $reactions_containing{$compound}->{$reaction} = 1;
   }
}
my %singletons_in_reaction;
foreach my $compound (keys(%reactions_containing))
{
    my $reactH = $reactions_containing{$compound};
    my @reactions = keys(%$reactH);
    if (@reactions == 1)
    {
        $singletons_in_reaction{$reactions[0]}->{$compound} = 1;
    }
}

foreach my $couplet (@couplets)
{
    my ($reaction,$row) = @$couplet;
    if (my $singletonsH = $singletons_in_reaction{$reaction})
    {
        my @singletons = sort keys(%$singletonsH);
        foreach my $compound (@singletons)
        {
            print join("\t",(@$row,$compound)),"\n";
        }
    }
}
