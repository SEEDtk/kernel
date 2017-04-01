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

=head1 Organize Shrub Quality Data

    organize_shrub_quality.pl [ options ]

This script takes output relating to Shrub genomes from the L<collate_checkm.pl> script and converts it to the
standard format produced by such scripts as L<package_report.pl>.

=head2 Parameters

There are no positional parameters. The standard input should be the output file from L<collate_checkm.pl>. The columns
of this file are

=over 4

=item 1

Genome ID

=item 2

Genome Name

=item 3

SciKit Coarse Score

=item 4

SciKit Fine Score

=item 5

CheckM Taxonomic Classification

=item 6

CheckM Completeness Score

=item 7

CheckM Contamination Score

=item 8

CheckM Heterogeneity Score

=back

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=over 4

=item core

Only include core genomes.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('',
        Shrub::script_options(),
        ScriptUtils::ih_options(),
        ['core', 'only include core genomes'],
       );
# Get the options.
my $coreOnly = $opt->core;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Loop through the input.
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    my ($genome, $name, $skCoarse, $skFine, $cmTaxon, $cmComplete, $cmContam) = split /\t/, $line;
    # Get the genome data from the shrub.
    my @contigData = $shrub->GetAll('Genome Contig', 'Genome(id) = ?', [$genome],
            'Genome(core) Genome(dna-size) Contig(length)');
    # Compute the genome source.
    my $core = $contigData[0][0];
    my $totLen = $contigData[0][1];
    my $source = ($core ? 'CoreSEED' : 'Shrub');
    # Only proceed if we found the genome and it matches our criteria.
    if ($totLen && (! $coreOnly || $core)) {
        # Get the genome metrics.
        my $totLen = $contigData[0][1];
        if (! $totLen) {
            print "ERROR in $genome.\n";
        }
        # Compute the contig-related data.
        my @lens = map { $_->[2] } @contigData;
        my $metrics = compute_metrics(\@lens, $totLen);
        # Output the data line.
        print join("\t", $source, $genome, $name, scalar(@contigData), $totLen, $metrics->{complete},
                $metrics->{N50}, $metrics->{N70}, $metrics->{N90}, $genome, $name, $skCoarse, $skFine,
                $cmComplete, $cmContam, $cmTaxon) . "\n";
    }
}


##Return a hash of metrics about this GTO. The metrics returned will include N50, N70, N90, total DNA length, and
## probable completeness.
sub compute_metrics {
    my ($lens, $totLen) = @_;
    # This will be the return hash.
    my %retVal;
    # Save the total length.
    $retVal{totlen} = $totLen;
    # Create a hash of threshholds.
    my %thresh = (N50 => 0.5 * $totLen, N70 => 0.7 * $totLen, N90 => 0.9 * $totLen);
    # Sort the contig lengths from longest to shortest.
    my @lens = sort { $b <=> $a } @$lens;
    # We accumulate the contig length as we go through the sorted list until we break a
    # threshold.
    my $cumul = 0;
    for my $len (@lens) {
        $cumul += $len;
        for my $type (keys %thresh) {
            if ($cumul >= $thresh{$type}) {
                # We have the desired metric. Save it in the return array.
                $retVal{$type} = $len;
                # Insure we don't test for it again.
                delete $thresh{$type};
            }
        }
    }
    # Check for completeness.
    $retVal{complete} = (($retVal{N70} >= 20000 && $totLen >= 300000) ? 1 : 0);
    # Return the hash.
    return \%retVal;
}
