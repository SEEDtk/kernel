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
use KmerDb;
use FastA;
use Stats;
use Shrub;

=head1 Create Signature DNA Kmers for Representative Genomes

    build_sample_sigs.pl [ options ] repDir

This script creates a signature kmer database for representative genomes. The genomes must all be stored in the Shrub
database.  The signature kmer database will be written to the file C<sigs.>I<NN>C<.json>,
where I<NN> is the kmer size.

=head2 Parameters

The positional parameter is the name of the representative-genome directory.

The command-line options are those in L<Shrub/script_options> to select the database, plus the following.

=over 4

=item kmer

Kmer size to use, in base pairs. The default is C<20>.

=item maxFound

The maximum number of times a kmer can be found before it is considered common. Common kmers are removed from the hash. The
default is C<10>.

=item indicative

If specified, the kmers will be indicative, rather than signature.  Normally, only kmers unique to a genome will be put in
the database. In indicative mode, all uncommon kmers are kept.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('repDir', Shrub::script_options(),
        ['kmer|K=i', 'kmer size to use', { default => 20 }],
        ['indicative|I', 'indicative instead of signature kmers'],
        ['maxFound|m=i', 'maximum number of occurrences before a kmer is considered common', { default => 10 }],
        );
my $stats = Stats->new();
# Get the representative genome directory.
my ($repDir) = @ARGV;
if (! $repDir) {
    die "No representative genome directory specified.";
} elsif (! -s "$repDir/6.1.1.20.dna.fasta") {
    die "$repDir does not have a DNA FASTA file.";
}
# Get the options.
my $K = $opt->kmer;
my $maxFound = $opt->maxfound;
my $indicMode = $opt->indicative;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Create the kmer database.
my $kmerDb = KmerDb->new(kmerSize => $K, maxFound => $maxFound);
# Read in the genome IDs. Compute the names and DNA files.
my %dnaFiles;
open(my $gh, "<$repDir/complete.genomes") || die "Could not open complete.genomes: $!";
while (! eof $gh) {
    my $line = <$gh>;
    if ($line =~ /^(\d+\.\d+)/) {
        my $genome = $1;
        my ($gData) = $shrub->GetAll('Genome', 'Genome(id) = ?', [$genome], 'contig-file name');
        if (! $gData) {
            die "$genome not found in database.";
        } else {
            $kmerDb->AddGroup($genome, $gData->[1]);
            $dnaFiles{$genome} = $gData->[0];
        }
        $stats->Add(genomeIn => 1);
    }
}
print "Genome names and DNA files stored.\n";
close $gh;
# Loop through the DNA files.
for my $genome (sort keys %dnaFiles) {
    my $name = $kmerDb->name($genome);
    my $dnaFile = "$FIG_Config::shrub_dna/$dnaFiles{$genome}";
    print "Processing $genome ($name) from $dnaFile.\n";
    $stats->Add(genomeScanned => 1);
    my $fh = FastA->new($dnaFile);
    while ($fh->next) {
        $kmerDb->AddSequence($fh->left);
        $stats->Add(sequenceIn => 1);
    }
}
# Now we finalize depending on the mode.
if ($indicMode) {
    print "Finalizing indicators.\n";
    $kmerDb->Finalize();
} else {
    print "Finalizing discriminators.\n";
    $kmerDb->ComputeDiscriminators();
}
# Write the output file.
my $outFile = "$repDir/sigs.$K.json";
print "Saving database as $outFile.\n";
$kmerDb->Save($outFile);
print "All done.\n" . $stats->Show();
