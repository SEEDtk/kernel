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
use Bin::ContigBlast;
use TetraMap;
use Bin::Analyze;
use Bin;
use Stats;
use Loader;
use SeedUtils;
use FastA;
use Data::Dump;

=head1 Create Bins From Community Contigs (Algorithm 5)

    bin2_generate.pl [ options ] sampleDir workDir

Process the crAss output for community samples to create bins from the contigs. We use the following algorithm.

=over 4

=item 1

Build bins from all the incoming contigs.

=item 2

Find locations of the default universal role C<PhenSyntTrnaAlph> in the contigs. Each hit will be to a specific
representative genome.

=item 3

Form a blast database of the universal roles in the representative genomes determined in step (2).

=item 4

Blast the contigs against the blast database formed in (3). The closest hit in each contig determines the bin.

=back

=head2 Parameters

There are two positional parameters

=over 4

=item 1

The name of a directory containing the sample. The sample's contigs must be in a file called C<contigs.fasta> and
the vector file in a file called C<output.contigs2reads.txt> in this directory. The directory may also contain
an optional list of excluded reference genomes-- C<exclude.tbl>.

=item 2

The name of the directory to contain the output and intermediate files. If this parameter is omitted, the input
directory is assumed.

=back

The command-line options are the following.

=over 4

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are used for computing
the identity of the bins. The default is C<400>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are used for
computing the identity of the bins. The default is C<4>.

=item binLenFilter

If specified, the minimum length of a contig that can be placed into a bin. The default is C<400>.

=item tetra

Tetranucleotide computation scheme. The default is C<dual>.

=over 8

=item raw

Each tetranucleotide and its reverse compliment is counted.

=item fancy

Each tetranucleotide and its reverse compliment is counted unless the reverse compliment is the same.

=item dual

Each tetranucleotide is counted once, and is considered identical to its reverse compliment.

=back

=item force

If C<all>, all of the data and intermediate files will be recomputed. If C<parms>, the initial bin files will be
reused, but everything else will be recomputed. (Use this last one if you have changed the parameters but not the
input data.) If omitted, nothing will be recomputed.

=item bindb

A directory containing protein FASTA files used for blasting. This directory is built by L<bin2_protein_database.pl>.

=item maxE

The maximum acceptable E-value. The default is C<1e-20>.

=item binMaxE

The maximum acceptable E-value for binning (as opposed to computation of reference genomes). The default is C<1e-20>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=item statistics-file

If specified, the name of a file into which the statistics should be saved.

=back

=head2 Input Files

The following files are expected in the sample input directory.

=over 4

=item contigs.fasta

This is a DNA FASTA file containing the contig ID and sequence for each contig to be binned.

=item output.contigs2reads.txt

This is a tab-delimited file. The first line is a header line and is discarded. Each subsequent line should contain
(0) a contig ID and (1) one or more coverage amounts (forming a vector of coverages by sample). Commonly, there is only
one coverage amount for every contig.

=back

=head2 Working Files

This script produces intermediate working files that are re-used if they already exist. Many files are output in the
working directory by L<Bin::Blast>.

The C<sample.fasta> file contains all of the sample sequences in FASTA form.

The C<reduced.fasta> file contains the sample sequences that pass the initial length and coverage filters.

The C<ref.genomes.scores.tbl> file contains the best hit for each seed contig along with its blast score.
It is a tab-delimited file containing (0) the seed contig ID, (1) the ID of the reference genome, (1) the percent identity
blast score, and (2) the scientific name of the reference genome.

The C<contigs.bins> file contains bin information for each input contig, in bin exchange format.
It does not include universal roles or reference genomes.

The C<proteins.ref.fa> file contains the universal proteins for the reference genomes. It will be used as a blast
database.

The C<bins.json> file contains the output bins, in json format.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir workDir',
                ['lenFilter=i',    'minimum contig length for seed protein search', { default => 400 }],
                ['covgFilter=f',   'minimum contig mean coverage for seed protein search', { default => 4}],
                ['tetra=s',        'tetranucleotide counting algorithm', { 'default' => 'dual' }],
                ['force=s',        'force re-creation of all intermediate files'],
                ['bindb=s',        'name of the directory containing the binning protein FASTA files',
                                    { default => "$FIG_Config::global/Bin5" }],
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-20 }],
                ['binMaxE=f',      'maximum acceptable e-value for binning blast hits', { default => 1e-20 }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
                ['statistics-file=s', 'save statistics data to this file'],
        );
# Turn off buffering for stdout.
$| = 1;
# Create the loader object.
my $loader = Loader->new();
# Get the statistics object.
my $stats = $loader->stats;
# Save the filter options.
my $lenFilter = $opt->lenfilter;
my $covgFilter = $opt->covgfilter;
# Check the file names.
my ($sampleDir, $workDir) = @ARGV;
if (! $sampleDir) {
    die "Sample directory name missing.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
my ($contigFile, $vectorFile) =
        map { "$sampleDir/$_" } qw(contigs.fasta output.contigs2reads.txt);
if (! -f $contigFile) {
    die "Contig file $contigFile not found.";
} elsif (! -f $vectorFile) {
    die "Vector file $vectorFile not found.";
}
# Check the working directory.
if (! $workDir) {
    $workDir = $sampleDir;
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.";
}
my ($sampleFastaFile, $reducedFastaFile, $binFile) = map { "$workDir/$_" } qw(sample.fasta reduced.fasta contigs.bins);
# Verify the bin data directory.
my $bindb = $opt->bindb;
if (! -d $bindb) {
    die "Binning protein directory $bindb not found.";
} elsif (! -s "$bindb/seedProtein.fa") {
    die "Binning protein directory $bindb missing seed protein file.";
}
# The $force flag will be set to TRUE if we want to force file creation at any point.
my $forceType = $opt->force // 'none';
my $force = ($forceType eq 'all');
# This hash will contain all the contigs, it maps each contig ID to a bin object describing the contig's properties.
my %contigs;
# Do we already have the initial contig bins?
if ($force || ! -s $reducedFastaFile || ! -s $binFile) {
    # We must process the raw contigs to create the contig bin objects and the sample FASTA file.
    # Create the tetranucleotide object.
    my $tetra = TetraMap->new($opt->tetra);
    # Now loop through the contig input file. We also save a copy of the contig sequences.
    print "Processing tetranucleotide data.\n";
    my $fh = $loader->OpenFasta(contig => $contigFile);
    open(my $ofh, '>', $sampleFastaFile) || die "Could not open FASTA output file: $!";
    my $fields = $loader->GetLine(contig => $fh);
    while (defined $fields) {
        my ($contig, undef, $seq) = @$fields;
        # Get the sequence length.
        my $len = length $seq;
        # Create a new bin for this contig.
        my $bin = Bin->new($contig);
        $bin->set_len($len);
        $stats->Add(contigLetters => $len);
        # Save it in the contig hash.
        $contigs{$contig} = $bin;
        # Save a copy of the sequence.
        print $ofh ">$contig\n$seq\n";
        # Compute the tetranucleotide vector.
        my $contigTetra = $tetra->ProcessString($seq);
        $bin->set_tetra($contigTetra);
        # Get the next line.
        $fields = $loader->GetLine(contig => $fh);
    }
    close $ofh;
    # Now all the bins are created. Next we get the coverage information. This is computed from the vector file.
    print "Reading coverage vector file.\n";
    my $vh = $loader->OpenFile(coverage => $vectorFile);
    # Throw away the first line. It contains headers.
    $fields = $loader->GetLine(coverage => $vh);
    # Loop through the rest of the file. There will be a mix of nodes and other lines. We only keep the nodes, which
    # have contig IDs matching what was in the contig FASTA.
    while (! eof $vh) {
        $fields = $loader->GetLine(coverage => $vh);
        my ($contigID, @coverages) = @$fields;
        # Get this contig's bin object.
        my $bin = $contigs{$contigID};
        if (! $bin) {
            $stats->Add(coverageLineSkipped => 1);
        } else {
            # Store the coverage vector in this contig.
            $bin->set_coverage(\@coverages);
            $stats->Add(coverageLineStored => 1);
        }
    }
    # Close the vector input file.
    close $vh;
    # Create a file of the input contigs and output the reduced FASTA file. This requires re-reading the
    # sample FASTA file.
    open(my $kfh, '>', $reducedFastaFile) || die "Could not open FASTA output file: $!";
    open(my $bfh, '>', $binFile) || die "Could not open contig bin output file: $!";
    $fh = $loader->OpenFasta(contig => $contigFile);
    $fields = $loader->GetLine(contig => $fh);
    while (defined $fields) {
        my ($contig, undef, $seq) = @$fields;
        my $bin = $contigs{$contig};
        if ($bin) {
            $bin->WriteContig($bfh);
            # Do we keep this contig for BLASTing against the seed protein?
            if ($bin->len < $lenFilter) {
                if ($bin->meanCoverage < $covgFilter) {
                    $stats->Add(contigRejectedBoth => 1);
                } else {
                    $stats->Add(contigRejectedLen => 1);
                }
            } elsif ($bin->meanCoverage < $covgFilter) {
                $stats->Add(contigRejectedCovg => 1);
            } else {
                # Yes. Save it.\
                print $kfh ">$contig\n$seq\n";
                $stats->Add(contigKept => 1);
            }
        }
        $fields = $loader->GetLine(contig => $fh);
    }
    # Force creation of all future files.
    $force = 1;
} else {
    # We can read the bin objects from the contigs.bin file and the FASTA file is already in place.
    my $binList = Bin::ReadContigs($binFile);
    for my $bin (@$binList) {
        $contigs{$bin->contig1} = $bin;
    }
}
# Turn on forcing if force = parms.
$force ||= ($forceType eq 'parms');
# Compute the blast parameters.
my $maxE = $opt->maxe;
my $minlen = $opt->minlen;
my $gap = $opt->gap;
my $binMaxE = $opt->binmaxe // $maxE;
# Get the name of the reference genome output file.
my $refGenomeFile = "$workDir/ref.genomes.scores.tbl";
# This will contain the tuples describing the reference genomes for each contig hit.
my %refTuples;
# Do we already have the reference genomes?
if ($force || ! -s $refGenomeFile) {
    # No, we must blast to find them.
    print "Searching for bins using seed protein.\n";
    my $blastMgr = Bin::ContigBlast->new("$bindb/seedProtein.fa", { maxE => $maxE, minLen => $minlen, tempDir => $workDir,
            gap => $gap, stats => $stats, logH => \*STDOUT });
    my $matches = $blastMgr->FindAllMatches($reducedFastaFile);
    # Open the output file.
    open(my $oh, '>', $refGenomeFile) || die "Could not open $refGenomeFile: $!";
    for my $contig (keys %$matches) {
        my $match = $matches->{$contig};
        # Write this match to the output file.
        my (undef, undef, undef, $score, $genome, $name) = @$match;
        print $oh join("\t", $contig, $genome, $score, $name) . "\n";
        # Save the reference tuple.
        $refTuples{$contig} = [$genome, $score, $name];
    }
} else {
    # Here we have our reference genomes already. Load the tuples from the file.
    print "Reading reference genomes from $refGenomeFile: $!";
    open(my $ih, '<', $refGenomeFile) || die "Could not open $refGenomeFile: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        my ($contig, $genome, $score, $name) = split /\t/, $line;
        $refTuples{$contig} = [$genome, $score, $name];
    }
}
# Now we begin the real binning. We are going to loop through the input FASTA file. For each contig that has a hit,
# we create a bin out of it for the specified reference genome.  All other contigs are written to a FASTA file for
# blasting.
my $unplacedFastaFile = "$workDir/unplaced.fasta";
open(my $fh, '>', $unplacedFastaFile) || die "Could not open $unplacedFastaFile: $!";
# This will track the bin ID for each reference genome. If more than one contig hits the same genome, the contigs
# are combined into a single bin.
my %gBin;
# This will map genome IDs to names.
my %gName;
# This will hold the result bins. Each will start with the bin object for a single contig.
my %bins;
# Loop through the contigs.
my $unCount = 0;
print "Collating contigs for binning.\n";
my $gh = FastA->new($sampleFastaFile);
while ($gh->next) {
    my $contig = $gh->id;
    my $tuple = $refTuples{$contig};
    # Is this contig already in a bin?
    if ($tuple) {
        # Yes. Get the genome ID.
        my $genome = $tuple->[0];
        # Is there already a bin for this genome?
        my $contig0 = $gBin{$genome};
        if ($contig0) {
            # Yes. Merge the bins.
            $bins{$contig0}->Merge($contigs{$contig});
            $stats->Add(binMerge => 1);
        } else {
            # No. Create a new bin for this genome.
            $gBin{$genome} = $contig;
            $bins{$contig} = $contigs{$contig};
            $gName{$genome} = $tuple->[2];
            $stats->Add(binCreated => 1);
            $bins{$contig}->add_ref($genome);
            # Compute the bin name.
            my ($genus, $species, $strain) = split /\s+/, $tuple->[2];
            if ($species eq 'sp.') {
                $species .= " $strain";
            }
            my $name = "$genus $species clonal population";
            my ($taxon) = split /\./, $genome;
            $bins{$contig}->set_name($name, $taxon);
            print "Bin created for contig $contig using $genome: $name.\n";
        }
    } else {
        # This contig is not in a bin. Save it to the FASTA file.
        $gh->Write($fh);
        $stats->Add(contigCopied => 1);
        $unCount++;
    }
}
undef $gh; close $fh;
# Now we create the universal protein database.
my $uniProtFile = "$workDir/proteins.ref.fa";
open(my $ph, '>', $uniProtFile) || die "Could not open $uniProtFile: $!";
for my $genome (sort keys %gBin) {
    print "Copying $genome.fa.\n";
    my $fh = FastA->new("$bindb/$genome.fa");
    while ($fh->next) {
        $fh->Write($ph, "$genome\t$gName{$genome}");
    }
}
close $ph;
# Now create the blast manager and blast the unbinned contigs against the universal proteins.
print "Blasting contigs against universal proteins.\n";
my $blastMgr = Bin::ContigBlast->new($uniProtFile, { maxE => $binMaxE, minLen => $minlen, tempDir => $workDir,
            gap => $gap, stats => $stats, logH => \*STDOUT });
my $matches = $blastMgr->FindAllMatches($unplacedFastaFile);
my $matchCount = scalar keys %$matches;
print "Matches found for $matchCount of $unCount contigs.\n";
# Loop through the matches. Each match will be for a contig that is not in a bin (that is, it is a unit), and
# the genome ID will identify a unique bin via the %gBin hash.
for my $contig (keys %$matches) {
    my $match = $matches->{$contig};
    my $genome = $match->[4];
    my $binContig = $gBin{$genome};
    $bins{$binContig}->Merge($contigs{$contig});
    $stats->Add(contigBinned => 1);
}
# Sort the bins and create the output file.
my @sorted = sort { Bin::cmp($a, $b) } values %bins;
print "Writing bins to output.\n";
open(my $oh, ">$workDir/bins.json") || die "Could not open bins.json file: $!";
for my $bin (@sorted) {
    $bin->Write($oh);
}
close $oh;
# Check to see if the user wants statistics saved.
if ($opt->statistics_file) {
    if (open(my $sh, ">", $opt->statistics_file)) {
        print $sh $stats->Show();
        close($sh);
    } else {
        warn "Cannot write statistics file " . $opt->statistics_file . ": $!";
    }
}
# All done.
print "All done.\n" . $stats->Show();
