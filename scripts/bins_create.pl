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
use TetraMap;
use Bin::Analyze;
use Bin;
use Stats;
use Loader;
use SeedUtils;

=head1 Create Bins From Community Contigs (Algorithm 3)

    bins_create.pl [ options ] sampleDir workDir

Process the crAss output for community samples to create bins from the contigs.

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
is C<1000>.

=item covgFilter

Minimum mean coverage amount for a contig to be considered meaningful. Only meaningful contigs are retained. The
default is C<10>.

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

=item unifile

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column,
occurrence counts in the second, and role names in the third.

=item force

If C<all>, all of the data and intermediate files will be recomputed. If C<parms>, the initial bin files will be
reused, but everything else will be recomputed. (Use this last one if you have changed the parameters but not the
input data.) If omitted, nothing will be recomputed.

=item seedrole

The ID of the universal role to use for seeding the bin assignment. The default is C<PhenTrnaSyntAlph>.

=item seedgenome

The ID of the genome to use for seeing the bin assignment. The default is E coli K12-- C<83333.1>. To
specify more than one genome, use a comma-delimited list (e.g. C<83333.1,224324.1,300852.3>).

=item maxE

The maximum acceptable E-value. The default is C<1e-20>.

=item refMaxE

The maximum acceptable E-value for blasting to determine the best reference genome for a seed contig. Each seed
contig eventually forms a bin. The default is C<1e-10>.

=item binMaxE

The maximum acceptable E-value for blasting to determine which bin a contig belongs in. The default is C<1e-30>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=item dna

Use DNA blasting to match reference genomes to sample contigs. This is the default.

=item prot

Use protein blasting to match reference genomes to sample contigs.

=item unassembled

Skip the final assembly step. The result will be that we have identified the bins, but we have not put contigs
into them.

=back

=head2 Working Files

This script produces intermediate working files that are re-used if they already exist. Many files are output in the
working directory by L<Bin::Blast>.

The C<sample.fasta> file contains the sequences for the meaningful contigs in the sample.

The C<rejected.fasta> file contains the sequences for the non-meaningful contigs in the sample.

The C<bins.found.tbl> file contains the locations hit by the universal protein used to seed the process.

The C<ref.genomes.tbl> file contains the reference genome for each bin.

The C<ref.genomes.scores.tbl> file contains the blast score for each reference genome and each bin.

The C<contigs.bins> file contains bin information for each contig considered to be meaningful, in bin exchange format.
It does not include universal roles or reference genomes.

The C<contigs.ref.bins> file contains bin information for each contig considered to be meaningful, in bin exchange format.
It includes reference genomes and (at the current time) a crude approximation of universal roles.

The C<bins.json> file contains the output bins, in json format.

The C<bins.report.txt> file contains a summary of the output bins.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir workDir',
                Shrub::script_options(),
                ['tetra=s',        'tetranucleotide counting algorithm', { 'default' => 'dual' }],
                ['unifile=s',      'universal role file', { default => "$FIG_Config::global/uni_roles.tbl" }],
                ['lenFilter=i',    'minimum contig length', { default => 1000 }],
                ['covgFilter=f',   'minimum contig mean coverage', { default => 10}],
                ['force=s',        'force re-creation of all intermediate files'],
                ['seedrole|R=s',   'ID of the universal role to seed the bins', { default => 'PhenTrnaSyntAlph' }],
                ['seedgenome|G=s', 'ID of the genome to seed the bins', { default => '83333.1,36870.1,224308.1,1148.1,64091.1,69014.3,83332.12,115711.7,187420.1,224326.1,243273.1,4932.3' }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-20 }],
                ['refMaxE=f',      'maximum acceptable e-value for reference genome blast hits', { default => 1e-10 }],
                ['binMaxE=f',      'maximum acceptable e-value for bin selection blast hits', { default => 1e-30 }],
                ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
                ['refsPerBin=i',   'number of reference genomes to keep per bin', { default => 1 }],
                ['unassembled',    'identify the bins, but do not assign contigs and create the final output file'],
                ['mode' => hidden => { one_of => [ ['dna|n' => 'use DNA blasting'], ['prot|p' => 'use protein blasting'] ] }]
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
my ($reducedFastaFile, $rejectedFastaFile, $binFile) = map { "$workDir/$_" } qw(sample.fasta rejected.fasta contigs.bins);
# Connect to the database.
print "Connecting to database.\n";
my $shrub = Shrub->new_for_script($opt);
# This will be set to TRUE if we want to force file creation at any point.
my $forceType = $opt->force // 'none';
my $force = ($forceType eq 'all');
# Get the seeding parameters.
my $prot = $opt->seedrole;
my $genome = $opt->seedgenome;
# Compute the mode.
my $blastMode = $opt->mode // 'dna';
if ($blastMode eq 'prot') {
    $blastMode = 'p';
} else {
    $blastMode = 'n'
}
# This hash will contain all the contigs, it maps each contig ID to a bin object describing the contig's properties.
my %contigs;
if ($force || ! -s $reducedFastaFile || ! -s $binFile) {
    # We must process the raw contigs to create the contig bin objects and the reduced FASTA file.
    # Create the tetranucleotide object.
    my $tetra = TetraMap->new($opt->tetra);
    # Now loop through the contig input file, filtering by length.
    print "Processing tetranucleotide data.\n";
    my $fh = $loader->OpenFasta(contig => $contigFile);
    my $fields = $loader->GetLine(contig => $fh);
    while (defined $fields) {
        my ($contig, undef, $seq) = @$fields;
        # Get the sequence length.
        my $len = length $seq;
        # Is this sequence long enough to be meaningful?
        if ($len < $lenFilter) {
            $stats->Add(contigTooShort => 1);
        } else {
            $stats->Add(contigLongEnough => 1);
            # Yes. Create a new bin for this contig.
            my $bin = Bin->new($contig);
            $bin->set_len($len);
            $stats->Add(contigLetters => $len);
            # Save it in the contig hash.
            $contigs{$contig} = $bin;
            # Compute the tetranucleotide vector.
            my $contigTetra = $tetra->ProcessString($seq);
            $bin->set_tetra($contigTetra);
        }
        # Get the next line.
        $fields = $loader->GetLine(contig => $fh);
    }
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
    # Discard low-coverage bins.
    for my $contig (keys %contigs) {
        my $bin = $contigs{$contig};
        if ($bin->meanCoverage < $covgFilter) {
            delete $contigs{$contig};
            $stats->Add(contigDeletedForCoverage => 1);
        }
    }
    print "Final count is " . scalar(keys %contigs) . " sample contigs kept.\n";
    # Create a file of the selected contigs. We need to re-open the input file and extract the
    # contigs we want to keep.
    open(my $ofh, '>', $reducedFastaFile) || die "Could not open FASTA output file: $!";
    open(my $xfh, '>', $rejectedFastaFile) || die "Could not open FASTA save file: $!";
    open(my $bfh, '>', $binFile) || die "Could not open contig bin output file: $!";
    $fh = $loader->OpenFasta(contig => $contigFile);
    $fields = $loader->GetLine(contig => $fh);
    while (defined $fields) {
        my ($contig, undef, $seq) = @$fields;
        if (exists $contigs{$contig}) {
            print $ofh ">$contig\n$seq\n";
            $contigs{$contig}->WriteContig($bfh);
        } else {
            print $xfh ">$contig\n$seq\n";
        }
        $fields = $loader->GetLine(contig => $fh);
    }
    # Force creation of all future files.
    $force = 1;
} else {
    # We can read the bin objects from the contigs.bin file and the FASTA files are already in place.
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
my $bMaxE = $opt->binmaxe // $maxE;
# Create the blaster.
my $blaster = Bin::Blast->new($shrub, $workDir, $reducedFastaFile, uniRoles => $opt->unifile,
        maxE => $maxE, minlen => $opt->minlen, gap => $opt->gap);
# First, we need the list of bins and the locations where they hit the primary universal protein.
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
    my $ih = $loader->OpenFile(binsFound => $binsFoundFile);
    while (my $binFields = $loader->GetLine(binsFound => $ih)) {
        my ($contig, $begin, $dir, $len) = @$binFields;
        $matches->{$contig} = BasicLocation->new($contig, $begin, $dir, $len);
    }
}
# Next, we need a hash of contig IDs to reference genomes. This is only for the contigs that represent the
# starting bins. We are essentially finding the best N reference genome DNA matches for each of the
# protein fragments identified in the first step.
my $contigHash;
# Do we have a reference genome file.
my $refGenomeFile = "$workDir/ref.genomes.tbl";
my $refScoreFile = "$workDir/ref.genomes.scores.tbl";
if ($force || ! -s $refGenomeFile) {
    # No. Create a hash mapping each contig ID to the DNA sequence representing the hit. We do this by reading
    # the reduced FASTA file and applying the matches hash.
    my $seqHash = $loader->GetDNA($matches, $reducedFastaFile);
    # Now BLAST against a database of the seed protein.
    $contigHash = $blaster->MatchProteins($seqHash, $prot, 1, $rMaxE);
    # Save the contig list to the reference-genomes files.
    open(my $oh, ">$refGenomeFile") || die "Could not open reference genome output file: $!";
    open(my $sh, ">$refScoreFile") || die "Could not open referecne genome score file: $!";
    for my $contig (sort keys %$contigHash) {
        my @genomes;
        for my $genomeData (@{$contigHash->{$contig}}) {
            my ($genome, $score) = @$genomeData;
            print $sh join("\t", $contig, $genome, $score) . "\n";
            push @genomes, $genome;
        }
        print $oh join("\t", $contig, @genomes) . "\n";
        # Remove the scores from the contig hash.
        $contigHash->{$contig} = \@genomes;
        $stats->Add('refGenomes-lineOut' => 1);
    }
    $force = 1;
} else {
    # Read the hash from the reference genome file.
    my $ih = $loader->OpenFile(refGenomes => $refGenomeFile);
    while (my $refFields = $loader->GetLine(refGenomes => $ih)) {
        my ($contig, @genomes) = @$refFields;
        $contigHash->{$contig} = \@genomes;
    }
}
# Our final list of bins goes in here.
my @binList;
# Create a bin for each reference genome found, and fill it with the starter contigs.
# We also put these bins in the bin list.
my %binHash;
for my $contig (keys %$contigHash) {
    my $rgList = $contigHash->{$contig};
    my ($rg) = @$rgList;
    my $bin = $contigs{$contig};
    $bin->add_ref($rg);
    my $gBin = $binHash{$rg};
    if ($gBin) {
        $gBin->Merge($bin);
        $stats->Add(starterBinCombined => 1);
    } else {
        $binHash{$rg} = $bin;
        push @binList, $bin;
        $stats->Add(starterBinCreated => 1);
    }
}
my @refGenomes = sort keys %binHash;
# The next step is to assign reference genomes to the bins. We cache this information in a file.
my $refBinFile = "$workDir/contigs.ref.bins";
if ($force || ! -s $refBinFile) {
    # Blast the reference genomes against the community contigs to assign them.
    my $subStats = $blaster->Process(\%contigs, \@refGenomes, type => $blastMode, maxE => $bMaxE);
    $stats->Accumulate($subStats);
    # Checkpoint our results.
    open(my $oh, ">$refBinFile") || die "Could not open augmented contig bin output file: $!";
    for my $contig (sort keys %contigs) {
        $contigs{$contig}->WriteContig($oh);
    }
    $force = 1;
} else {
    # Read the augmented contig information.
    %contigs = ();
    my $binList = Bin::ReadContigs($refBinFile);
    for my $bin (@$binList) {
        $contigs{$bin->contig1} = $bin;
    }
}
# Now we want to run through the augmented contigs, forming them into bins by reference genome.
if ($opt->unassembled) {
    print "Skipping assembly step.\n";
} else {
    print "Assembling bins.\n";
    for my $bin (values %contigs) {
        my ($refGenome) = $bin->refGenomes;
        if (! $refGenome) {
            push @binList, $bin;
            $stats->Add(contigNoRefGenome => 1);
        } elsif (exists $binHash{$refGenome}) {
            my $oldBin = $binHash{$refGenome};
            if ($oldBin ne $bin) {
                $oldBin->Merge($bin);
                $stats->Add(contigMerged => 1);
            }
        } else {
            $binHash{$refGenome} = $bin;
            push @binList, $bin;
            $stats->Add(contigNewBin => 1);
        }
    }
    # Sort the bins and create the initial report.
    my @sorted = sort { Bin::cmp($a, $b) } @binList;
    print "Writing bins to output.\n";
    open(my $oh, ">$workDir/bins.json") || die "Could not open bins.json file: $!";
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
    my $analyzer = Bin::Analyze->new(totUnis => $totUnis, minUnis => (0.8 * $totUnis));
    $analyzer->BinReport($rh, $uniRoles, \@sorted);
    close $rh;
}
# All done.
print "All done.\n" . $stats->Show();
