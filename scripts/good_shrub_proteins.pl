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
use GPUtils;
use File::Copy::Recursive;
use Stats;
use Shrub;

=head1 Produce Report on Seed Protein in Shrub Genomes

    good_proteins.pl [ options ] outDir

This script produces two files-- an amino-acid FASTA file of all occurrences of a selected protein in the specified
genome packages, and a tab-delimited file listing each genome ID and its name.

=head2 Output Files

Two output files will be put in the designated output directory.

=over 4

=item complete.genomes

A tab-delimited file of genomes, each record consisting of (0) a genome ID and (1) a genome name.

=item 6.1.1.20.fasta

A FASTA file containing the amino acid sequences of the features containing the desired protein role. The
sequence ID will be the feature ID. The comment will be the genome ID.

=back

=head2 Parameters

The positional parameter is the name of a directory to contain the output files.

The standard input should contain a list of genomes to examine. The input should be tab-delimited, with the
genome ID in the first column.

The command-line options are those in L<Shrub/script_options> plus those in L<ScriptUtils/ih_options> for specifying
the standard input, and the following.

=over 4

=item protein

The ID of the desired protein role. The default is C<PhenTrnaSyntAlph>.

=item minlen

The minimum acceptable length for the protein. The default is 209.

=item maxlen

The maximum acceptable length for the protein. The default is 485.

=item col

The (1-based) input column containing the genome IDs. Use C<0> for the last column. The default is C<1>.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(), ScriptUtils::ih_options(),
        ['protein=s', 'protein role description', { default => 'PhenTrnaSyntAlph'}],
        ['minlen=i', 'minimum protein length', { default => 209 }],
        ['maxlen=i', 'maximum protein length', { default => 485 }],
        ['col=i', 'genome ID column (1-based)', { default => 1 }],
        );
# Get the positional parameters.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) ||
        die "Could not create $outDir: $!";
}
# Get the options.
my $role = $opt->protein;
my $min = $opt->minlen;
my $max = $opt->maxlen;
my $col = $opt->col - 1;
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the standard input.
my $ih = ScriptUtils::IH($opt->input);
# Open the output files.
open(my $fh, '>', "$outDir/6.1.1.20.fasta") || die "Could not open FASTA output file: $!";
open(my $oh, '>', "$outDir/complete.genomes") || die "Could not open genome output file: $!";
# Get all of the genomes in the input.
print "Scanning input file.\n";
my %gHash;
while (! eof $ih) {
    my $line = <$ih>;
    $line =~ s/[\r\n]+$//;
    $stats->Add(lineIn => 1);
    my @cols = split /\t/, $line;
    my $genome = $cols[$col];
    my ($name) = $shrub->GetFlat('Genome', 'Genome(id) = ?', [$genome], 'name');
    if ($name) {
        $gHash{$genome} = $name;
        $stats->Add(genomeFound => 1);
    }
}
print scalar(keys %gHash) . " genomes found.\n";
# Get all of the seed proteins in the Shrub.
print "Reading Shrub seed proteins.\n";
my %prots;
my $q = $shrub->Get('Role2Function Function2Feature Feature2Protein Protein',
        'Role2Function(from-link) = ? AND Function2Feature(security) = ?', [$role, 0],
        'Feature2Protein(from-link) Protein(sequence)');
while (my $pData = $q->Fetch()) {
    my ($fid, $seq) = $pData->Values('Feature2Protein(from-link) Protein(sequence)');
    if ($fid =~ /fig\|(\d+\.\d+)/) {
        my $genome = $1;
        if ($gHash{$genome}) {
            push @{$prots{$genome}}, [$fid, $seq];
            $stats->Add(protFound => 1);
        }
    }
}
# Loop through the incoming genomes.
for my $genome (sort keys %gHash) {
    print "Searching $genome for seed proteins:  ";
    $stats->Add(genomesIn => 1);
    my $protList = $prots{$genome};
    if (! $protList) {
        $stats->Add(genomesNoSeed => 1);
        print "none found.\n";
    } elsif (scalar(@$protList) > 1) {
        my $found = scalar(@$protList);
        print "$found found. GENOME SKIPPED.\n";
        $stats->Add(("genomes$found" . "Seed") => 1);
    } else {
        print "1 found.";
        $stats->Add(genomes1Seed => 1);
        # Get the feature ID.
        my $fid = $protList->[0][0];
        # Check the protein length.
        my $aa = $protList->[0][1];
        my $aaLen = length($aa);
        if ($aaLen < $min) {
            print " Length $aaLen < $min. GENOME SKIPPED.\n";
            $stats->Add(protTooShort => 1);
        } elsif ($aaLen > $max) {
            print " Length $aaLen > $max. GENOME SKIPPED.\n";
            $stats->Add(protTooLong => 1);
        } else {
            print "\n";
            # We found the protein. Output the genome.
            print $oh "$genome\t$gHash{$genome}\n";
            $stats->Add(genomesOut => 1);
            # Output the protein.
            print $fh ">$fid $genome\n$aa\n";
            $stats->Add(proteinsOut => 1);
        }
    }
}
print "All done.\n" . $stats->Show();