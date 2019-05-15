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
use P3Utils;
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

=item gLabel

If specified, the FASTA sequence labels will be genome IDs instead of feature IDs, and the comments will be
feature IDs and genome names.

=item gList

If specified, the name of a tab-delimited file containing genome IDs in the first column. Only genome IDs
from the list will be used.

=back

=cut

# Insure output is up-to-date.
$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('dbName',
                Shrub::script_options(),
                ['seedrole|R=s',   'ID of the universal role to seed the bins', { default => 'PhenTrnaSyntAlph' }],
                ['prot|p',         'output amino acid sequences instead of DNA sequences'],
                ['gLabel|G',       'use genome labeling'],
                ['gList=s',        'specify allowed genomes']
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
# Get the inclusion list.
my $inclusions;
if ($opt->glist) {
    my $file = $opt->glist;
    open(my $gh, '<', $opt->glist) || die "Could not open genome list: $!";
    while (! eof $gh) {
        my $line = <$gh>;
        if ($line =~ /^(\d+.\d+)/) {
            $inclusions->{$1} = 1;
        }
    }
    print scalar(keys %$inclusions) . " genomes in inclusion list.\n";
}
# Determine the labeling scheme.
my $glabel = $opt->glabel;
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
    if ($inclusions && ! $inclusions->{$genomeID}) {
        $stats->Add(genomeNotIncluded => 1);
    } elsif (! $taxID) {
        $stats->Add(genomeNoTaxID => 1);
    } else {
        $stats->Add(genomeGood => 1);
        $genomes{$genomeID} = $genomeName;
    }
}
# Now we have a list of the genomes we can use.
print scalar(keys %genomes) . " genomes with good taxonomy information.\n";
# Release the genome list memory.
@genomes = ();
# Now we need to get all the features for the seed role.
my ($roleData) = $shrub->GetAll("Role", 'Role(id) = ?', [$opt->seedrole], 'Role(description) Role(checksum)');
my ($roleName, $roleCheck) = @$roleData;
# Clean the role name for PATRIC and ask for features.
$roleName =~ s/\([^)]+\)/ /g;
$roleName =~ s/\s+/ /g;
$roleName =~ s/^\s+//;
$roleName =~ s/\s+$//;
print "Role name is $roleName.\n";
my $prots = P3Utils::get_data($p3, feature => [['eq', 'product', $roleName]], ["patric_id", "genome_id", "product", $seqType]);
# This will hold our protein batch.
# Loop through the proteins, keeping the ones for the good genomes.
my $totCount = scalar(@$prots);
print "$totCount proteins found.\n";
my $protCount = 0;
for my $prot (@$prots) {
    my ($protID, $genomeID, $product, $seq) = @$prot;
    my $genomeName = $genomes{$genomeID};
    my $protCheck = RoleParse::Checksum($product);
    $stats->Add(protFound => 1);
    if (! $genomeName) {
        # Relevant genome is ambiguously classified.
        $stats->Add(protBadGenome => 1);
    } elsif (! $protID) {
        # Not a PATRIC feature.
        $stats->Add(protNotPatric => 1);
    } elsif ($protCheck ne $roleCheck) {
        # Wrong role returned.
        $stats->Add(protWrongRole => 1);
    } else {
        # Compute the label/comment combination for this protein.
        my $header = ($glabel ? "$genomeID $protID" : "$protID $genomeID" ) . "\t$genomeName";
        # Process the protein.
        if (! $seq) {
            $stats->Add(missingSequence => 1);
        } else {
            print $oh ">$header\n$seq\n";
            $stats->Add(protOut => 1);
        }
    }
    $protCount++;
    print "$protCount of $totCount proteins processed.\n" if $protCount % 10000 == 0;
}
print "All done:\n" . $stats->Show();

