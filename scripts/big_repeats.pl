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
use ScriptUtils;
use BlastUtils;

=head1 Find Repeat Regions in Genomes

    big_repeats.pl [ options ] > repeats.txt

Find regions that appear to be big repeats at the DNA level. This can be done by looking for multiple copies of
identical DNA within a signle genome or looking for instances of large repeats maintained as a blast database.

=head2 Parameters

The contigs to be examined for repeat regions are specified via the standard input (or the C<--input> option) in
FASTA format.

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item genome

ID of a genome in the Shrub database to be searched for repeat regions. If this option is specified, the
standard input ignored.

=item minIden

The minimum percent identity value for a region to be considered a repeat. The default is C<95>, indicating
a 95% match.

=item minLen

The minimum length for a region to be considered a repeat. The default is C<100>.

=item blastDB

If this option is specified, a repeat is defined as a similarity against an entry in the specified
blast database.

=back

=head2 Output Format

The output is to the standard output, and is tab-delimited with 4 columns. Each record represents a single
repeat.

=over 4

=item 0

region length

=item 1

percent identity

=item 2

location string for the first region

=item 3

location string for the repeat region

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('parms', Shrub::script_options(), ScriptUtils::ih_options(),
        ['genome|g=s',    'id of a genome to search'],
        ['minIden|m=i',   'minimum percent identity', { default => 95 }],
        ['minLen|l=i',    'minimum region length', { default => 100 }],
        ['blastDB|b=i',   'repeat region blast database']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Now we need the query file.. This goes in the following variable.
my $queryFile;
# Do we have a genome ID?
if ($opt->genome) {
    # Yes. Get the genome's FASTA file.
    $queryFile = $shrub->genome_fasta($opt->genome);
} else {
    # No genome ID. Do we have an input file?
    if ($opt->input) {
        # Yes. Return it.
        $queryFile = $opt->input;
    } else {
        # No. Make a query file out of the standard input.
        $queryFile = BlastUtils::get_query(\*STDIN);
        if (! $queryFile) {
            die "Could not create FASTA query file from standard input.";
        }
    }
}
# Next we need the database. This goes in the following variable.
my $blastDB;
my $haveDB;
# Do we have a user-specified blast database?
if ($opt->blastdb) {
    # Yes, so use it.
    $blastDB = $opt->blastdb;
    $haveDB = 1;
} else {
    # No. Make a blast database out of the query file.
    $blastDB = BlastUtils::get_db($queryFile, 'blastn', $FIG_Config::temp);
    if (! $blastDB) {
        die "Could not create BLAST database out of $queryFile.";
    }
}
# This hash will track matches we have already seen.
my %seen;
# Now we run the BLAST. If we find a match longer than the minimum length that does not overlap itself,
# we treat it as a repeat region and write it to the output. Note that we pre-filter on the length and
# percent identity.
my $matches = BlastUtils::blast($queryFile, $blastDB, 'blastn',
                { maxE => 1e-20, minIden => $opt->miniden / 100, minlen => $opt->minlen,
                  outForm => 'hsp', includeSelf => 1 });
for my $match (@$matches) {
    # Get the match locations. We use this to insure the match is unique and is non-overlapping.
    my ($loc1, $loc2) = ($match->qloc, $match->sloc);
    # Only proceed if this is not a match of a sequence against itself.
    if ($haveDB || ! $loc1->OverlapLoc($loc2) ) {
        # Sort the match parameters.
        if (BasicLocation::Cmp($loc1, $loc2) < 0) {
            ($loc1, $loc2) = ($loc2, $loc1);
        }
        # Form a hash string to insure this match is unique.
        my $string1 = $loc1->String;
        my $string2 = $loc2->String;
        my $hashString = join("\t", $string1, $string2);
        # Only proceed if it is not.
        if (! $seen{$hashString}) {
            # Insure we don't see it again.
            $seen{$hashString} = 1;
            # Write it out.
            print join("\t", $match->n_mat, $match->pct, $string1, $string2) . "\n";
        }
    }
}
