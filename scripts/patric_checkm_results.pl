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
use Shrub;
use SeedUtils;
use Stats;

=head1 Process PATRIC CheckM Data

    patric_checkm_results.pl [ options ] contigData checkmDataDir

This script creates a skeleton package quality report from a directory of CheckM output files and a tab-delimited file of
contig data. The SciKit data will need to be filled in later.

The standard columns for a package quality report are

=over 4

=item 1

Sample name (always C<PATRIC>)

=item 2

Genome ID

=item 3

Genome scientific name

=item 4

Number of contigs

=item 5

Number of DNA base pairs

=item 6

C<1> if the genome is mostly complete, else C<0>

=item 7

The contig N50 score.

=item 8

The contig N70 score.

=item 9

The contig N90 score.

=item 10

ID of the genome (same as column 2)

=item 11

Name of the genome (same as column 3)

=item 12

SciKit coarse evaluation score (always 0)

=item 13

SciKit fine evaluation score (always 0)

=item 14

CheckM evaluation completeness (percent)

=item 15

CheckM evaluation contamination (percent)

=item 16

CheckM taxonomy classification

=back

The scientific name is taken by querying the database using the taxonomy ID from the contig file. The data on DNA sizes and
completeness is computed from the contig lengths.

=head2 Parameters

The positional parameters are the name of the contig input file and the name of a directory containing CheckM output files.

The contig input file is a tab-delimited file with a header line. The columns are (0) a genome ID, (1) the taxonomy ID for
the genome, (2) a contig name, and (3) a contig length. The file must be sorted by genome ID.

The CheckM files should all have names in the form I<XXXXXX.X>C<.checkm.out> where I<XXXXXX.X> is a genome ID.

The command-line options are those found in L<Shrub/script_options>.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('contigData checkmDataDir', Shrub::script_options(),
        );
# Initialize a statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Verify the parameters.
my ($contigData, $checkmDataDir) = @ARGV;
if (! $contigData) {
    die "No contig input file specified.";
} elsif (! -s $contigData) {
    die "Contig data file $contigData missing or empty.";
} elsif (! $checkmDataDir) {
    die "No checkM data directory specified.";
} elsif (! -d $checkmDataDir) {
    die "CheckM data direcory $checkmDataDir missing or invalid.";
}
# This hash will map taxonomy IDs to taxonomic grouping names.
my %taxon;
# This will contain the current genome and taxon IDs.
my ($genomeID, $taxonID) = ('', '');
# This will contain the list of contig lengths.
my @lens;
# Read the input file.
open(my $ih, "<$contigData") || die "Could not open contig input file: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($genome, $taxon, $contig, $len) = split /\t/, $line;
    $stats->Add(contigLineIn => 1);
    # Do we have a new genome?
    if ($genome ne $genomeID) {
        # Yes. Insure the old one (if any) is processed.
        if ($genomeID) {
            OutputGenome($genomeID, $taxonID, \@lens);
        }
        # Start the new genome.
        ($genomeID, $taxonID) = ($genome, $taxon);
        @lens = ($len);
    } else {
        # Old genome. Add this length.
        push @lens, $len;
    }
}
# Process the residual genome.
OutputGenome($genomeID, $taxonID, \@lens);
# All done.
print STDERR "Statistics for this run:\n" . $stats->Show();

## Output the data for a genome.
sub OutputGenome {
    my ($genomeID, $taxonID, $lens) = @_;
    $stats->Add(genomeProcessed => 1);
    # We will store the checkM data in these variables.
    my ($complete, $contam, $taxon);
    # Look for the checkM data file.
    if (! open(my $ih, "<$checkmDataDir/$genomeID.checkm.out")) {
        print STDERR "Could not access checkM data for $genomeID: $!\n";
        $stats->Add(noCheckmData => 1);
    } else {
        # Find the checkM data in the file.
        while (! eof $ih && ! defined $complete) {
            my $line = <$ih>;
            if ($line =~ /\s+\d+\.\d+\s+(\S+)/) {
                $taxon = $1;
                my @parts = split /\s\s+/, $line;
                ($complete, $contam) = @parts[12,13];
            }
        }
        if (! defined $complete) {
            print STDERR "Invalid checkM file format for $genomeID.\n";
            $stats->Add(badCheckmData => 1);
        }
    }
    # Compute the taxonomic classification.
    my ($taxName) = $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(id) = ?', [$taxonID], 'scientific-name');
    if (! $taxName) {
        print STDERR "No taxonomic data found for $genomeID using $taxonID.\n";
        $stats->Add(badTaxonID => 1);
    }
    if (defined $complete && $taxName) {
        # Here we have all the external data we need. Compute the metrics.
        my $metricsH = SeedUtils::compute_metrics($lens);
        print join("\t", 'PATRIC', $genomeID, $taxName, scalar(@$lens), $metricsH->{totlen},
                $metricsH->{complete}, $metricsH->{N50}, $metricsH->{N70}, $metricsH->{N90},
                $genomeID, $taxName, 0, 0, $complete, $contam, $taxon) . "\n";
        $stats->Add(genomeOut => 1);
    }
}