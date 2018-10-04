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
use Stats;
use File::Copy::Recursive;

=head1 Create FASTA Files for Protein-Based Binning (Algorithm 5)

    bin2_protein_database.pl [ options ] outDir

This script creates the FASTA files for the latest (generation 5) binning algorithm. A file is created for
the seed protein (PhenTrnaSyntAlph) in all of the Shrub genomes, and then one is created for all the universal
proteins in each genome. Not absolutely all of the genomes in Shrub are used: only the well-behaved core genomes
and the prokaryotic PATRIC genomes are included. Optionally, a list of genomes to be used can be selected.

The comment for each sequence will be the genome ID followed by a tab and the genome name.

=head2 Parameters

The positional parameter is the name of the directory in which the output files are to be created. There will be
approximately 4000 files, which is a lot.

The command-line options are those found in L<Shrub/script_options>, for selecting the Shrub database, plus the
following:

=over 4

=item clear

If specified, the output directory will be erased if it already exists. If not specified, an existing output directory
will cause an error.

=item gList

If specified, a list of the genomes to be used. The file must be tab-delimited, with the genome ID in the first
column and the genome name in the second. This overrides the criteria for selecting genomes described above.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ['clear', 'erase output directory if it exists'],
        ['gList=s', 'tab-delimited file of genome IDs and names to use']
        );
my $stats = Stats->new();
# Check the parameters.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (-f $outDir) {
    die "Invalid output directory $outDir specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} elsif ($opt->clear) {
    print "Erasing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase $outDir: $!";
} else {
    die "$outDir already exists. Use --clear to erase it.";
}
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
print "Connected to Shrub.\n";
# First, get a hash table of the genomes of interest.
my %genomes;
my $gFile = $opt->glist;
if (! $gFile) {
    print "Loading genomes from database.\n";
    %genomes = map { $_->[0] => $_->[1] } $shrub->GetAll("Genome",
            'Genome(prokaryotic) = ? AND (Genome(core) = ? OR Genome(well-behaved) = ?)',
            ['1','0','1'], [qw(id name)]);
} else {
    print "Reading genomes from $gFile.\n";
    open(my $gh, "<$gFile") || die "Could not open genome file: $!";
    while (! eof $gh) {
        my $line = <$gh>;
        if ($line =~ /^(\d+\.\d+)\t(.+)/) {
            $genomes{$1} = $2;
        }
    }
}
my $gCount = scalar keys %genomes;
print "$gCount genomes selected for processing.\n";
# Get all of the seed proteins.
my @prots = $shrub->GetAll("Function2Feature Feature2Protein Protein",
        'Function2Feature(from-link) = ? AND Function2Feature(security) = ?',
        ['PhenTrnaSyntAlph','1'],
        [qw(Feature2Protein(from-link) Protein(sequence))]);
print scalar(@prots) . " seed proteins found.\n";
# Open the output FASTA file.
open(my $fh, ">$outDir/seedProtein.fa") || die "Could not open output FASTA file: $!";
for my $prot (@prots) {
    my ($fid, $seq) = @$prot;
    my $genome = SeedUtils::genome_of($fid);
    my $name = $genomes{$genome};
    if (! $name) {
        $stats->Add(seedProtBadGenome => 1);
    } else {
        print $fh ">$fid $genome\t$name\n$seq\n";
        $stats->Add(seedProtOut => 1);
    }
}
close $fh; undef $fh;
print "seedProtein.fa created.\n";
my $count = 0;
# Now we create a FASTA file for all the universal proteins in each genome.
for my $genome (sort keys %genomes) {
    my $name = $genomes{$genome};
    $count++;
    print "Processing $genome ($count of $gCount): $name.\n";
    $stats->Add(genomeProcessed => 1);
    open($fh, ">$outDir/$genome.fa") || die "Could not open FASTA file for $genome: $!";
    # Get all the universal proteins for this genome.
    @prots = $shrub->GetAll("Genome2Feature Feature Feature2Function Function AND Feature Protein",
            'Genome2Feature(from-link) = ? AND Feature2Function(security) = ? AND Function(universal) = ?',
            [$genome,'1','1'], [qw(Feature(id) Protein(sequence))]);
    # Loop through them, producing output.
    for my $prot (@prots) {
        my ($fid, $seq) = @$prot;
        print $fh ">$fid $genome\t$name\n$seq\n";
        $stats->Add(protProcessed => 1);
    }
    close $fh; undef $fh;
}
print "All done.\n" . $stats->Show();
