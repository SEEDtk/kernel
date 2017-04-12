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
use P3DataAPI;
use Stats;

=head1 Generate Seed Protein Database for Binning

    bins_protein_database.pl [ options ] dbName

This script generates the seed protein database for binning. The seed protein database is a BLAST database
that contains each instance of the seed protein found in PATRIC. It is used to find which PATRIC genome is
closest to each seed protein instance, and that in turn determines the reference genome for a bin.

The database is a FASTA file. The ID is a FIG feature ID. The comment contains two tab-delimited fields-- a SEED-style genome ID
and a genome name. Only genomes whose genome name matches the name of their associated genus are included, because only these are
considered unambiguously identified.

=head2 Parameters

The sole positional parameter is the name of the protein database. This must be the full name of a FASTA file to be
generated; however, any BLAST database work files with the same name will be deleted before this program starts.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item seedrole

The ID of the universal role to use for seeding the bin assignment. The default is C<PhenTrnaSyntAlph>.

=item prot

If specified, the FASTA will contain amino acid sequences instead of DNA sequences.

=back

=cut

# Insure output is up-to-date.
$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dbName',
                Shrub::script_options(),
                ['seedrole|R=s',   'ID of the universal role to seed the bins', { default => 'PhenTrnaSyntAlph' }],
                ['prot',           'output amino acid sequences instead of DNA sequences'],
        );
# Get the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Prepare the database.
my ($dbName) = @ARGV;
if (! $dbName) {
    die "No output database name specified.";
} elsif (-f $dbName) {
    print "Deleting old copy of $dbName.\n";
    my @files = glob("$dbName.*");
    push @files, $dbName;
    for my $file (@files) {
        unlink $file;
        $stats->Add(fileDeleted => 1);
    }
}
# Determine the output sequence type.
my $seqType = 'na_sequence';
if ($opt->prot) {
    $seqType = 'aa_sequence';
}
# Open the FASTA file for output.
open(my $oh, ">$dbName") || die "Could not open output FASTA file: $!";
# Enable access to PATRIC from Argonne.
$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME} = 0;
# Connect to PATRIC.
my $p3 = P3DataAPI->new();
my @genomes = $p3->query("genome", ["ne", "genome_id", 0],
        ["select", "genome_id", "genome_name", "taxon_id"]);
print scalar(@genomes) . " genomes found in PATRIC.\n";
# Loop through the PATRIC genomes, saving the good ones.
my %genomes;
for my $genome (@genomes) {
    $stats->Add(genomes => 1);
    # Get the genome's data.
    my $taxID = $genome->{taxon_id};
    my $genomeID = $genome->{genome_id};
    my $genomeName = $genome->{genome_name};
    if (! $taxID) {
        $stats->Add(genomeNoTaxID => 1);
    } else {
        # Extract the genus from the genome name.
        my ($genus) = split ' ', $genomeName;
        # Get the genus name from the taxonomy ID.
        my $taxName;
        while (! $taxName && $taxID ne 1) {
            my ($taxData) = $shrub->GetAll('TaxonomicGrouping IsInTaxonomicGroup', 'TaxonomicGrouping(id) = ?', [$taxID],
                'TaxonomicGrouping(scientific-name) TaxonomicGrouping(type) IsInTaxonomicGroup(to-link)');
            if (! $taxData) {
                # Here we have an invalid taxonomy ID.
                $stats->Add(genomeNoTaxId => 1);
                $taxID = 1;
            } else {
                my ($newName, $type, $nextID) = @$taxData;
                $stats->Add(genusSearch => 1);
                if ($type eq 'genus') {
                    $stats->Add(genusFound => 1);
                    $taxName = $newName;
                } else {
                    $taxID = $nextID;
                }
            }
        }
        if (! $taxName) {
            $stats->Add(genomeNoGenus => 1);
        } elsif ($taxName ne $genus) {
            $stats->Add(genomeBadGenus => 1);
        } else {
            $stats->Add(genomeGood => 1);
            $genomes{$genomeID} = $genomeName;
        }
    }
}
# Now we have a list of the genomes we can use.
print scalar(keys %genomes) . " genomes with good genus information.\n";
# Release the genome list memory.
@genomes = ();
# Now we need to get all the features for the seed role.
my ($funcName) = $shrub->GetFlat("Function", 'Function(id) = ?', [$opt->seedrole], 'Function(description)');
print "Function name is $funcName.\n";
# Clean the function name for PATRIC and ask for features.
$funcName =~ s/[()]/ /g;
my @prots = $p3->query("genome_feature", ["select", "patric_id", "genome_id"],
        ["eq", "annotation", "PATRIC"], ["eq", "product", qq("$funcName")]);
# This will hold our protein batch.
my %prots;
my $batchCount = 0;
my $protCount = 0;
# Loop through the proteins, keeping the ones for the good genomes.
print scalar(@prots) . " proteins found.\n";
for my $prot (@prots) {
    my $genomeID = $prot->{genome_id};
    my $protID = $prot->{patric_id};
    my $genomeName = $genomes{$genomeID};
    $stats->Add(protFound => 1);
    if (! $genomeName) {
        # Relevant genome is ambiguously classified.
        $stats->Add(protBadGenome => 1);
    } else {
        $prots{$protID} = "$genomeID\t$genomeName";
        $protCount++;
        $batchCount++;
        if ($batchCount >= 100) {
            ProcessProteins($oh, $p3, $stats, \%prots, $seqType);
            $batchCount = 0;
            %prots = ();
            print "$protCount proteins processed.\n";
        }
    }
}
if ($batchCount > 0) {
    ProcessProteins($oh, $p3, $stats, \%prots, $seqType);
}
print "All done:\n" . $stats->Show();

sub ProcessProteins {
    my ($oh, $p3, $stats, $prots, $seqType) = @_;
    my @protList = $p3->query("genome_feature", ["select", "patric_id", $seqType],
            ["in", "patric_id", '(' . join(",", keys %$prots) . ')']);
    for my $prot (@protList) {
        my $label = $prot->{patric_id};
        my $comment = $prots->{$label};
        my $seq = $prot->{$seqType};
        if (! $seq) {
            $stats->Add(missingSequence => 1);
        } else {
            print $oh ">$label $comment\n$seq\n";
            $stats->Add(protOut => 1);
        }
    }
}