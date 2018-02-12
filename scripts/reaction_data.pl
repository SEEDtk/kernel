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
use Shrub;
use File::Copy::Recursive;
use ScriptUtils;

=head1 Produce Data on Chemical Reactions

    reaction_data.pl [ options ] outDir

This script dumps out two files detailing the reactions stored in the L<Shrub> database.

=head2 Output Files

The following files will be put into the output directory.

=over 4

=item reactions.txt

Each reaction is described by five lines in a tab-delimite file.

=over 8

=item 1

(0) the reaction ID and (1) C<Y> for reversible reactions and a null string otherwise.

=item 2

A tab-separated list of substrate compound IDs.

=item 3

A tab-separated list of product compound IDs.

=item 4

A printable version of the reaction.

=item 5

A constant delimiter C<//>.

=back

=item compounds.tbl

A tab-delimited file, one line per compound, containing (0) the compound ID and (1) the compound name.

=back

=head2 Parameters

The positional parameter is the name of the output directory for the files. If it does not exist, it will be created.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item genome

Specifies a genome ID. If specified, only reactions relevant to that genome will be output.

=back

=cut

use constant CONNECTORS => { '<' => '<=', '=' => '<=>', '>' => '=>' };

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ['genome=s', 'ID of genome to which reactions should be restricted']
        );
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (-f $outDir) {
    die "Invalid output directory $outDir.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive($outDir);
}
# Connect to the database.
print "Connecting to the database.\n";
my $shrub = Shrub->new_for_script($opt);
# Get the list of relevant reactions.
my %reactions;
if (! $opt->genome) {
    print "Getting all reactions.\n";
    %reactions = map { $_ => 1 } $shrub->GetFlat('Reaction', '', [], 'id');
} else {
    my $genome = $opt->genome;
    print "Getting reactions for $genome.\n";
    %reactions = map { $_ => 1 } $shrub->GetFlat('Genome2Feature Feature2Function Function2Role Role2Complex Complex2Reaction',
            'Genome2Feature(from-link) = ? AND Feature2Function(security) = ?', [$genome, 0], 'Complex2Reaction(to-link)');
}
print scalar(keys %reactions) . " reactions found.\n";
# The compounds will be kept in here. It maps IDs to names.
my %compounds;
# Open the main output file.
open(my $oh, ">$outDir/reactions.txt") || die "Could not open reactions.txt: $!";
my $count = 0;
# Loop through the reactions.
for my $rxn (sort keys %reactions) {
    my @formulaData = $shrub->GetAll('Reaction Reaction2Compound Compound', 'Reaction(id) = ?', [$rxn],
            "Reaction(direction) Reaction2Compound(product) Reaction2Compound(stoichiometry) Compound(id) Compound(label)");
    # Only proceed if we found something.
    if (@formulaData) {
        # This will contain the compounds. 0 is substrates, 1 is products.
        my @compounds = ({ }, { });
        # We accumulate the left and right sides of the formula separately in here.
        my @side = ([], []);
        # This will hold the connector.
        my $dir;
        for my $formulaDatum (@formulaData) {
            # Get the information for this compound.
            my ($direction, $product, $stoich, $cID, $cForm) = @$formulaDatum;
            # Save the compound for the compound file.
            $compounds{$cID} = $cForm;
            # Create the formula element.
            my $compound = ($stoich > 1 ? "$stoich*" : '') . $cForm;
            # Save it in the appropriate  formula queue.
            push @{$side[$product]}, $compound;
            # Save the compound ID for the product/substrate report.
            $compounds[$product]{$cID} = 1;
            # Save the direction.
            $dir //= CONNECTORS->{$direction};
        }
        # Construct the formula.
        my $formula = join(" $dir ", map { join(" + ", @$_) } @side);
        # Compute the reversible flag.
        my $revFlag = ($dir eq '<=>' ? 'Y' : '');
        # Write the reaction.
        print $oh join("\t", $rxn, $revFlag) . "\n";
        for my $key (0, 1) {
            print $oh join("\t", sort keys %{$compounds[$key]}) . "\n";
        }
        print $oh "$formula\n//\n";
        $count++; print "$count reactions processed.\n" if $count % 100 == 0;
    }
}
# Write out the compound report.
print "Writing compound report.\n";
close $oh; undef $oh;
open($oh, ">$outDir/compounds.tbl") || die "Could not open compounds.tbl: $!";
for my $cID (sort keys %compounds) {
    print $oh "$cID\t$compounds{$cID}\n";
}
print "All done.\n";