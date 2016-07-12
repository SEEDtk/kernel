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
use Bin::Blast;
use Bin::Tetra;
use Bin::Analyze;
use Bin::Kmer;
use Bin;
use Stats;
use Loader;
use SeedUtils;
use P3DataAPI;
use GenomeTypeObject;
use KmerDb;

=head1 Create Bins From Community Contigs (Algorithm 4)

    bins_generate.pl [ options ] sampleDir workDir

Process the crAss output for community samples to create bins from the contigs. We use the following algorithm.

=over 4

=item 1

Build bins from all the incoming contigs.

=item 2

Find locations of the default universal role C<PhenSyntTrnaAlph>, taken from a small set of sample genomes
in SEEDtk, in contigs with a certain minimum coverage (default 4) and length (default 400).

=item 3

Process the blast hits from step (2). BLAST them against a database of all the
PATRIC instances of the default universal roles. For each one, choose the best hit.
This is the bin's representative genome.

=item 4

Merge bins whose representative genomes belong to the same genus. Each bin will now have a
set of representative genomes associated with it and one or more initial contigs.

=item 5

Pull the representative genomes from PATRIC.  Build a protein kmer database (default length 12) from them. Each input
contig that contains more than a specified number of minimum number of occurrences of kmers from the same bin's
representative genomes (default 10 kmers) goes in that bin.

=item 6

Build a DNA kmer database (default length 50) from the binned contigs. Each unbinned contig with a matching kmer goes
in that bin.

=back

The kmer databases built in steps (5) and (6) only include kmers that appear in a single bin.

Unlike previous binning algorithms, this one does not output an estimate of universal roles in each bin.
One of the analysis scripts must be used to compute these.

=head2 Parameters

There are two positional parameters

=over 4

=item 1

The name of a directory containing the sample. The sample's contigs must be in a file called C<contigs.fasta> and
the vector file in a file called C<output.contigs2reads.txt> in this directory.

=item 2

The name of the directory to contain the output and intermediate files. If this parameter is omitted, the input
directory is assumed.

=back

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are retained. The default
is C<400>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are retained. The
default is C<4>.

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

=item seedrole

The ID of the universal role to use for seeding the bin assignment. The default is C<PhenTrnaSyntAlph>.

=item seedgenome

The ID of the genome to use for seeing the bin assignment. To specify more than one genome, use a
comma-delimited list. The default is C<83333.1,36870.1,224308.1,1148.1,64091.1,69014.3,83332.12,115711.7,187420.1,224326.1,243273.1,4932.3>.

=item seedfasta

The name of the BLAST database for the seed protein in the various PATRIC genomes. The default is
C<PhenTrnaSyntAlph.fa> in the global data directory.

=item maxE

The maximum acceptable E-value. The default is C<1e-20>.

=item refMaxE

The maximum acceptable E-value for blasting to determine the best reference genome for a seed contig. Each seed
contig eventually forms a bin. The default is C<1e-10>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=item unassembled

Skip the final assembly step. The result will be that we have identified the bins, but we have not put contigs
into them.

=item kmer

Kmer length for protein matching. The default is C<12>.

=item binstrength

The number of kmer matches required to place a contig into a bin. The default is C<10>.

=item danglen

Kmer length for placing unbinned contigs. The default is C<50>.

=item species

If specified, reference genomes will be grouped by genus and species instead of genus.

=back

=head2 Working Files

This script produces intermediate working files that are re-used if they already exist. Many files are output in the
working directory by L<Bin::Blast>.

The C<sample.fasta> file contains all of the sample sequences in FASTA form. 

The C<bins.found.tbl> file contains the locations hit by the universal protein used to seed the process.

The C<ref.genomes.scores.tbl> file contains the best hit for each seed contig along with its blast score.
It is a tab-delimited file containing (0) the seed contig ID, (1) the ID of the reference genome, (1) the percent identity
blast score, and (2) the scientific name of the reference genome.

The C<contigs.bins> file contains bin information for each input contig, in bin exchange format.
It does not include universal roles or reference genomes.

The C<kmers.json> file contains the bin-placement kmer database, in json format.

The reference genomes will be stored in L<GenomeTypeObject> json format in files named C<XXXXXXX.json>,
where I<XXXXXXX> is the genome ID.

The C<bins.kmer.json> file contains the first set of output bins (based on protein kmers) in json format.

The C<bins.unplaced.bin> file contains the unplaced contig bins from the protein kmer processing, in json format.

The C<unplaced.fa> file contains the contigs (in FASTA format) that were not placed by the protein kmer processing.

The C<bins.json> file contains the output bins, in json format.

The C<bins.report.txt> file contains a summary of the output bins.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir workDir',
                Shrub::script_options(),
                ['lenFilter=i',    'minimum contig length', { default => 400 }],
                ['covgFilter=f',   'minimum contig mean coverage', { default => 4}],
                ['tetra=s',        'tetranucleotide counting algorithm', { 'default' => 'dual' }],
                ['force=s',        'force re-creation of all intermediate files'],
                ['seedrole|R=s',   'ID of the universal role to seed the bins', { default => 'PhenTrnaSyntAlph' }],
                ['seedgenome|G=s', 'ID of the genome to seed the bins', { default => '83333.1,36870.1,224308.1,1148.1,64091.1,69014.3,83332.12,115711.7,187420.1,224326.1,243273.1,4932.3' }],
                ['seedfasta|F=s',  'BLAST database (or FASTA file) of seed protein in all genomes', { default => "$FIG_Config::global/PhenTrnaSyntAlph.fa"}],
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-20 }],
                ['refMaxE=f',      'maximum acceptable e-value for reference genome blast hits', { default => 1e-10 }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
                ['unassembled',    'identify the bins, but do not assign contigs and create the final output file'],
                ['kmer|k=i',       'kmer length for protein matches during binning', { default => 12 }],
                ['binstrength=i',  'number of kmer matches required to bin a contig', { default => 10 }],
                ['danglen=i',      'kmer length for unbinned-contig DNA matches', { default => 50 }],
                ['species',          'group by species instead of genus'],
        );
# Enable access to PATRIC from Argonne.
$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME} = 0;
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
# Verify the Seed Protein FASTA file.
my $seedFastaFile = $opt->seedfasta;
if (! -s $seedFastaFile) {
    die "Seed FASTA file $seedFastaFile not found.";
}
# Connect to the database.
print "Connecting to database.\n";
my $shrub = Shrub->new_for_script($opt);
# This will be set to TRUE if we want to force file creation at any point.
my $forceType = $opt->force // 'none';
my $force = ($forceType eq 'all');
# Get the seeding parameters.
my $prot = $opt->seedrole;
my $genome = $opt->seedgenome;
# This hash will contain all the contigs, it maps each contig ID to a bin object describing the contig's properties.
my %contigs;
# Do we already have the initial contig bins?
if ($force || ! -s $reducedFastaFile || ! -s $binFile) {
    # We must process the raw contigs to create the contig bin objects and the sample FASTA file.
    # Create the tetranucleotide object.
    my $tetra = Bin::Tetra->new($opt->tetra);
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
my $rMaxE = $opt->refmaxe // $maxE;
# Create the blaster.
my $blaster = Bin::Blast->new($shrub, $workDir, $reducedFastaFile, 
        maxE => $maxE, minlen => $opt->minlen, gap => $opt->gap);
# First, we need the list of bins and the locations where they hit the seed protein.
my $binsFoundFile = "$workDir/bins.found.tbl";
my $matches = {};
# Do we have a bins-found file?
if ($force || ! -s $binsFoundFile) {
    # No. we must search for the specified universal protein to create the initial bins.
    # Get the list of genomes.
    my @genomes = split /,/, $genome;
    my $string = $genome;
    if (scalar(@genomes) == 2) {
        $string = join(" and ", @genomes);
    } elsif (scalar(@genomes) > 2) {
        my @sorted = sort @genomes;
        my $last = pop @sorted;
        $string = join(", ", @sorted) . ", and $last";
    }
    print "Seeding bin process with $prot from $string.\n";
    $matches = $blaster->FindProtein(\@genomes, $prot);
    # Save the matches to a work file.
    open(my $oh, ">$binsFoundFile") || die "Could not open bins found output file: $!";
    for my $contig (sort keys %$matches) {
        my $match = $matches->{$contig};
        print $oh join("\t", $match->Contig, $match->Begin, $match->Dir, $match->Length) . "\n";
        $stats->Add('binsFound-lineOut' => 1);
    }
    $force = 1;
} else {
    # Yes. Read the initial bin information from the bins-found file.
    print "Reading bin definitions from $binsFoundFile.\n";
    my $ih = $loader->OpenFile(binsFound => $binsFoundFile);
    while (my $binFields = $loader->GetLine(binsFound => $ih)) {
        my ($contig, $begin, $dir, $len) = @$binFields;
        $matches->{$contig} = BasicLocation->new($contig, $begin, $dir, $len);
    }
}
# Next, we need a hash of contig IDs to reference genomes. This is only for the contigs that represent the
# starting bins. We find the best match for each region found in the previous step.
my $contigHash;
# This hash will contain the IDs and names of the reference genomes.
my %genomes;
# Do we have a reference genome file?
my $refScoreFile = "$workDir/ref.genomes.scores.tbl";
if ($force || ! -s $refScoreFile) {
    # No. Create a hash mapping each contig ID to the DNA sequence representing the hit. We do this by reading
    # the sample FASTA file and applying the matches hash.
    my $seqHash = $loader->GetDNA($matches, $reducedFastaFile);
    # Now BLAST against a database of the seed protein.
    $contigHash = $blaster->MatchProteins($seqHash, $prot, 1, $rMaxE, db => $seedFastaFile, type => 'dna');
    # Save the contig list to the reference-genome score file.
    open(my $sh, ">$refScoreFile") || die "Could not open reference genome score file: $!";
    for my $contig (sort keys %$contigHash) {
        for my $genomeData (@{$contigHash->{$contig}}) {
            my ($genome, $score, $name) = @$genomeData;
            print $sh join("\t", $contig, $genome, $score, $name) . "\n";
            $genomes{$genome} = $name;
        }
        $stats->Add('refGenomes-lineOut' => 1);
    }
    $force = 1;
} else {
    # Read the hash from the reference genome file.
    print "Reading reference genomes from $refScoreFile.\n";
    my $ih = $loader->OpenFile(refGenomes => $refScoreFile);
    while (my $refFields = $loader->GetLine(refGenomes => $ih)) {
        my ($contig, $genome, $score, $name) = @$refFields;
        push @{$contigHash->{$contig}}, [$genome, $score, $name];
        $genomes{$genome} = $name;
    }
}
# Our final list of bins goes in here.
my @binList;
# And this is the full list of reference genomes.
my %rg;
# Create a bin for each reference genome found, and fill it with the starter contigs.
# We also compute the best genome for the bin. This is the reference genome with the highest
# blast score.
my %binHash;
my %binBest;
for my $contig (keys %$contigHash) {
    # MatchProteins returns a list of results for each contig, but we have set up the parms to insure
    # each list is a singleton.
    my $genomeData = $contigHash->{$contig}[0];
    my ($genome, $score, $name) = @$genomeData;
    # Compute the title for this genome, depending on whether we are sorting on genus or genus/species.
    my ($genus, $species) = split ' ', $name;
    my $title = ($opt->species ? join(' ', $genus, $species) : $genus);
    my $bin = $contigs{$contig};
    $bin->add_ref($genome);
    $rg{$genome} = 1;
    my $gBin = $binHash{$title};
    if ($gBin) {
        $gBin->Merge($bin);
        $stats->Add(starterBinCombined => 1);
        my $oldScore = $binBest{$title}[1];
        if ($score > $oldScore) {
            $binBest{$title} = [$genome, $score];
        }
    } else {
        $binHash{$title} = $bin;
        push @binList, $bin;
        $stats->Add(starterBinCreated => 1);
        $binBest{$title} = [$genome, $score];
    }
}
my @refGenomes = sort keys %rg;
# Now we need GTOs for all these reference genomes.
my $p3 = P3DataAPI->new();
my @gtos;
# This hash tracks the GTOs that need to be saved to disk.
my %gtosToSave;
for my $refGenome (@refGenomes) {
    my $refGenomeFile = "$workDir/$refGenome.json";
    my $gto;
    if (! -f $refGenomeFile || $force) {
        print "Reading $refGenome from web API.\n";
        $gto = $p3->gto_of($refGenome);
        $gtosToSave{$refGenome} = $gto;
        $stats->Add(refGenomeFromWeb => 1);
        $force = 1;
    } else {
        print "Reading $refGenome from $refGenomeFile.\n";
        $gto = GenomeTypeObject->create_from_file($refGenomeFile);
        $stats->Add(refGenomeFromFile => 1);
    }
    if (! $gto) {
        die "$refGenome not found.";
    }
    push @gtos, $gto;
    $rg{$refGenome} = $gto;
}
# Get the name and taxonomy data for each bin.
for my $title (keys %binBest) {
    my $genome = $binBest{$title}[0];
    my $gto = $rg{$genome};
    $binHash{$title}->set_name("$title clonal population", $gto->{ncbi_taxonomy_id});
}
my $kmerDB;
# The next step is to assign contigs to the bins.
if ($opt->unassembled) {
    print "Skipping assembly step.\n";
    # Save the GTOs.
    SaveGTOs(\@refGenomes, \%gtosToSave, \@gtos, $workDir);
} else {
    print "Assembling bins.\n";
    # Check for the kmer database.
    my $kmerFile = "$workDir/kmers.json";
    if (! -s $kmerFile || $force) {
        # Create a new, blank kmer database.
        $kmerDB = Bin::Kmer->new(kmerSize => $opt->kmer, binStrength => $opt->binstrength);
        # Loop through the bins.
        for my $bin (@binList) {
            my $binID = $bin->contig1;
            print "Processing kmers for bin $binID.\n";
            $stats->Add(kmerBin => 1);
            # Loop through the genomes for this bin.
            my @genomes = $bin->refGenomes;
            for my $genome (@genomes) {
                my $gto = $rg{$genome};
                print "Computing kmers for $genome.\n";
                $stats->Add(kmerGenome => 1);
                $kmerDB->AddGenome($gto, $binID);
            }
        }
        # Compute the discriminating kmers.
        $kmerDB->Finalize();
        # Save the kmer database.
        print "Saving kmer database to $kmerFile.\n";
        $kmerDB->Save($kmerFile);
        # Denote the environment has changed.
        $force = 1;
    } else {
        # Read the kmer database.
        print "Reading kmer database from $kmerFile.\n";
        $kmerDB = Bin::Kmer->new(json => $kmerFile);
    }
    # Save the GTOs.
    SaveGTOs(\@refGenomes, \%gtosToSave, \@gtos, $workDir);
    # This will be a hash of the main bins.
    my %binHash;
    # Do we have a kmer-based version of the bins?
    my $kmerBinFile = "$workDir/bins.kmer.json";
    my $contigFile = "$workDir/bins.unplaced.bin";
    if (! $force && -s $kmerBinFile && -f $contigFile) {
        # Yes. Read it in.
        print "Reading bin status from $kmerBinFile\n";
        my $binsTemp = Bin::ReadBins($kmerBinFile);
        @binList = @$binsTemp;
        %binHash = map { $_->contig1 => $_ } @binList;
        %contigs = ();
        print "Reading contig status from $contigFile\n";
        $binsTemp = Bin::ReadContigs($contigFile);
        for my $bin (@$binsTemp) {
            $contigs{$bin->contig1} = $bin;
        }
    } else {
        # No. Clear the placed contigs from the contig hash.
        for my $bin (@binList) {
            my @contigs = $bin->contigs;
            $binHash{$bin->contig1} = $bin;
            for my $contig (@contigs) {
                delete $contigs{$contig};
            }
        }
        # Loop through the contigs, placing the unplaced ones.
        my ($contigCount, $placeCount) = (0, 0);
        print "Processing sample contigs.\n";
        my $fh = $loader->OpenFasta(sample => $sampleFastaFile);
        my $fields = $loader->GetLine(contig => $fh);
        while (defined $fields) {
            my ($contig, $comment, $sequence) = @$fields;
            if ($contigs{$contig}) {
                $contigCount++;
                my $binID = $kmerDB->ComputeBin($sequence, translate => 11);
                if ($binID) {
                    # We found a bin, so place this contig.
                    $stats->Add(placedContigKmer => 1);
                    my $bin = $binHash{$binID};
                    $bin->Merge($contigs{$contig});
                    delete $contigs{$contig};
                    $placeCount++;
                } else {
                    # No definite bin. Leave the contig unplaced.
                    $stats->Add(unplacedContigKmer => 1);
                }
                if ($contigCount % 1000 == 0) {
                    print "$contigCount contigs processed, $placeCount placed.\n";
                }
            }
            $fields = $loader->GetLine(contig => $fh);
        }
        # Checkpoint the bins.
        print "Writing unplaced contigs.\n";
        open(my $ch, '>', $contigFile) || die "Could not open contigs checkpoint file: $!";
        for my $contig (sort keys %contigs) {
            $contigs{$contig}->WriteContig($ch);
        }
        close $ch;
        print "Writing kmer-based bins.\n";
        my @sorted = sort { Bin::cmp($a, $b) } @binList;
        open(my $oh, ">$kmerBinFile") || die "Could not open bins.json checkpoint file: $!";
        for my $bin (@sorted) {
            $bin->Write($oh);
        }
        close $oh;
    }
    # Release the memory for the kmers.
    undef $kmerDB;
    # Now we need to look for mobile elements. We want to create a kmer database for long DNA kmers in the
    # binned contigs and use it to classify the unplaced contigs. While we are doing this, we put
    # the unplaced contigs in a FASTA file so we can process them faster. To facilitate this process,
    # we start with a hash that maps the IDs of the placed contigs to their bins.
    print "Cataloging binned contigs.\n";
    my %contig2Bin;
    for my $bin (@binList) {
        my $binID = $bin->contig1;
        my @contigs = $bin->contigs;
        for my $contig (@contigs) {
            $contig2Bin{$contig} = $binID;
        }
    }
    print "Building mobile element kmer database.\n";
    $kmerDB = Bin::Kmer->new(kmerSize => $opt->danglen, binStrength => 1);
    my $fh = $loader->OpenFasta(sample => $sampleFastaFile);
    my $unplacedFastaFile = "$workDir/unplaced.fasta";
    open(my $oh, '>', $unplacedFastaFile) || die "Could not open unplaced sequence output file: $!";
    my $fields = $loader->GetLine(contig2 => $fh);
    my ($contigCount, $placeCount) = (0, 0);
    while (defined $fields) {
        my ($contig, undef, $sequence) = @$fields;
        my $binID = $contig2Bin{$contig};
        if (! $binID) {
            print $oh ">$contig\n$sequence\n";
        } else {
            $kmerDB->AddSequence($binID, $sequence);
            $placeCount++;
        }
        $contigCount++;
        if ($contigCount % 1000 == 0) {
            print "$placeCount of $contigCount sequences processed for mobile element kmers.\n";
        }
        $fields = $loader->GetLine(contig2 => $fh);
    }
    print "Finalizing mobile element kmers from $placeCount sequences.\n";
    $kmerDB->Finalize();
    close $oh;
    # Now read the unplaced contigs and try to place them with the mobile elements.
    print "Processing unplaced contigs for mobile elements.\n";
    ($contigCount, $placeCount) = (0, 0);
    $fh = $loader->OpenFasta(unplaced => $unplacedFastaFile);
    $fields = $loader->GetLine(unplaced => $fh);
    while (defined $fields) {
        my ($contig, $comment, $sequence) = @$fields;
        $contigCount++;
        my $binID = $kmerDB->ComputeBin($sequence);
        if ($binID) {
            # We found a bin, so place this contig.
            $stats->Add(placedContigMobile => 1);
            my $bin = $binHash{$binID};
            $bin->Merge($contigs{$contig});
            $placeCount++;
        } else {
            # No definite bin. Leave the contig unplaced.
            $stats->Add(unplacedContigMobile => 1);
        }
        if ($contigCount % 1000 == 0) {
            print "$contigCount contigs processed, $placeCount placed.\n";
        }
        $fields = $loader->GetLine(unplaced => $fh);
    }
    # Sort the bins and create the initial report.
    my @sorted = sort { Bin::cmp($a, $b) } @binList;
    print "Writing bins to output.\n";
    undef $oh;
    open($oh, ">$workDir/bins.json") || die "Could not open bins.json file: $!";
    for my $bin (@sorted) {
        $bin->Write($oh);
    }
    close $oh;
    # Get the universal role information.
    print "Producing report.\n";
    my $uniRoles = $blaster->uni_hash;
    my $totUnis = $blaster->uni_total;
    # Produce the bin report.
    open(my $rh, ">$workDir/bins.report.txt") || die "Could not open report file: $!";
    my $analyzer = Bin::Analyze->new($shrub, totUnis => $totUnis, minUnis => (0.8 * $totUnis));
    $analyzer->SetGenomes(\%genomes);
    $analyzer->BinReport($rh, $uniRoles, \@sorted);
    close $rh;
}
# All done.
print "All done.\n" . $stats->Show();

# Save the GTOs to disk and release them from memory.
sub SaveGTOs {
    my ($refGenomes, $gtosToSave, $gtos, $workDir) = @_;
    for my $refGenome (@$refGenomes) {
        my $gto = $gtosToSave->{$refGenome};
        if ($gto) {
            my $gtoFileName = "$workDir/$refGenome.json";
            print "Saving genome to $gtoFileName.\n";
            $gto->destroy_to_file($gtoFileName);
        }
    }
    # Release the GTO memory.
    $gtosToSave = {};
    $gtos = [];
}
