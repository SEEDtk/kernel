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
use Bin::Kmers;
use Bin;
use Stats;
use Loader;
use SeedUtils;

=head1 Analyze Community Contigs

    bins_init.pl [ options ] sampleDir outputFile

Process the crAss output for community samples to create the L<bin exchange format|Bin/Bin Exchange Format> for the contigs.

=head2 Parameters

There are two positional parameters.

=over 4

=item 1

The name of a directory containing the sample. The sample's contigs must be in a file called C<contigs.fasta> and
the vector file in a file called C<output.contigs2reads.txt> in this directory.

=item 2

The name of the output file to contain the bin exchange format data for the contigs.

=back

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item minsim

Minimum permissible percent identity score for a BLAST result. The default is C<0.9>.

=item minlen

Minimum permissible length for a BLAST match. The default is C<150>.

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

=item workDir

Working directory to contain the intermediate files.

=item refs

The name of a file containing the IDs of the reference genomes. If omitted, all prokaryotic genomes will be
used. The file is tab-delimited, with the genome IDs in the first column.

=item unis

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column.

=item reps

The name of a file containing the IDs of the representative roles. If omitted, the file C<rep_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column.

=item force

If specified, the BLAST and KMER databases will be rebuilt even if they already exists.

=item minhits

The minimum nubmer of kmer hits required for a representative genome to be considered for BLASTing.

=item maxFound

The maximum number of kmer occurrences that can be found before a kmer is considered to common to be useful.

=item kmersize

The protein kmer size.

=back

=head2 Working Files

This script produces intermediate working files that are re-used if they already exist. Many files are output in the
working directory by L<Bin::Blast> and L<Bin::Kmers>. These are independent of the contigs input, so the same working
directory can be used for multiple communities.

The C<ref.counts.tbl> file contains the IDs of the reference genomes to be used for the analysis of the sample.
If it exists and the C<force> option is not specified, it is read to get the final reference genome IDs. Otherwise,
it is computed. This is stored in the sample directory.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir outputFile',
                Shrub::script_options(),
                ['minsim=f',     'minimum percent identity score for BLAST results', { 'default' => 0.9 }],
                ['minlen=i',     'minimum length for an acceptable BLAST match', { 'default' => 150 }],
                ['tetra=s',      'tetranucleotide counting algorithm', { 'default' => 'dual' }],
                ['workDir|D=s',  'working directory for intermediate files', { required => 1 }],
                ['refs=s',       'reference genome list file'],
                ['unis=s',       'universal role file', { default => "$FIG_Config::data/Global/uni_roles.tbl" }],
                ['reps=s',       'representative role file', { default => "$FIG_Config::data/Global/rep_roles.tbl" }],
                ['force',        'force a rebuild of the BLAST and KMER databases'],
                ['minhits=i',    'minimum number of kmer hits for a genome to be considered useful', { default => 400 }],
                ['maxfound=i',   'maximum number of kmer occurrences before it is considered common', { default => 10 }],
                ['kmersize|k=i', 'kmer size (protein)', { default => 10 }],
        );
# Create the loader object.
my $loader = Loader->new();
# Get the statistics object.
my $stats = $loader->stats;
# Save the force-flag.
my $force = $opt->force // 0;
# Check the working directory.
my $workDir = $opt->workdir;
if (! -d $workDir) {
    die "Invalid working directory $workDir.";
}
# Check the file names.
my ($sampleDir, $outputFile) = @ARGV;
if (! $sampleDir) {
    die "Sample directory name missing.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
my ($contigFile, $vectorFile) = map { "$sampleDir/$_" } qw(contigs.fasta output.contigs2reads.txt);
if (! -f $contigFile) {
    die "Contig file $contigFile not found.";
} elsif (! -f $vectorFile) {
    die "Vector file $vectorFile not found.";
} elsif (! $outputFile) {
    die "Output file name missing.";
}
# Connect to the database.
print "Connecting to database.\n";
my $shrub = Shrub->new_for_script($opt);
# Compute the reference genomes.
my $refs;
if ($opt->refs) {
    # Here we have a list file.
    $refs = $loader->GetNamesFromFile(inputGenomes => $opt->refs);
    print "Reference genomes read from file.\n";
} else {
    print "Computing reference genomes from database.\n";
    $refs = [ $shrub->GetFlat('Genome', 'Genome(prokaryotic) = ?', [1], 'id') ];
}
my $refCount = scalar @$refs;
print "$refCount reference genomes found.\n";
# This hash will contain all the contigs, it maps each contig ID to a bin object describing the contig's properties.
my %contigs;
# Create the tetranucleotide object.
my $tetra = Bin::Tetra->new($opt->tetra);
# Get the list of universal role IDs.
my $uniRoles = $loader->GetNamesFromFile(uniRoles => $opt->unis);
my $uniCount = scalar @$uniRoles;
print "$uniCount universal roles found.\n";
# Create the BLAST database.
print "Creating BLAST database.\n";
my $blaster = Bin::Blast->new($shrub, $workDir, $refs, $uniRoles, minsim => $opt->minsim, minlen => $opt->minlen, force => $force);
# Merge the statistics.
$stats->Accumulate($blaster->stats);
# Do we need to find the reference genomes?
my @realRefs;
my $refsFile = "$sampleDir/ref.counts.tbl";
if (-f $refsFile && ! $force) {
    # No, just read them in.
    my $refsList = $loader->GetNamesFromFile(finalRefGenomes => $refsFile);
    @realRefs = @$refsList;
} else {
    # Get the list of representative role IDs.
    my $repRoles = $loader->GetNamesFromFile(repRoles => $opt->reps);
    my $repCount = scalar @$repRoles;
    print "$repCount representative roles found.\n";
    # Create the KMER database.
    print "Creating KMER database.\n";
    my $kmers = Bin::Kmers->new($shrub, $workDir, $repRoles, force => $force, kmerSize => $opt->kmersize, minHits => $opt->minhits,
            maxFound => $opt->maxfound);
    # Compute the list of genomes to BLAST. First we find the relevant genomes.
    print "Computing reference genomes.\n";
    my $refGenomeH = $kmers->FindGenomes($contigFile);
    $stats->Add(closeRefGenome => scalar(keys %$refGenomeH));
    # Prepare to write what we decide to keep.
    open(my $oh, ">", $refsFile) || die "Could not open reference genome save file: $!";
    # Filter by the known representatives.
    for my $refGenome (@$refs) {
        my $count = $refGenomeH->{$refGenome};
        if ($count) {
            push @realRefs, $refGenome;
            $stats->Add(keptRefGenome => 1);
            print $oh "$refGenome\t$count\n";
        }
    }
}
# Now loop through the contig input file.
my $fh = $loader->OpenFasta(contig => $contigFile);
my $fields = $loader->GetLine(contig => $fh);
while (defined $fields) {
    my ($contig, undef, $seq) = @$fields;
    print "Processing contig $contig.\n";
    # Create a new bin for this contig.
    my $bin = Bin->new($contig);
    my $len = length $seq;
    $bin->set_len($len);
    $stats->Add(contigLetters => $len);
    # Save it in the contig hash.
    $contigs{$contig} = $bin;
    # Compute the tetranucleotide vector.
    my $contigTetra = $tetra->ProcessString($seq);
    $bin->set_tetra($contigTetra);
    # Compute the BLAST data.
    my ($genomesL, $rolesL) = $blaster->Process($seq, \@realRefs);
    $bin->add_ref(@$genomesL);
    if (! @$genomesL) {
        $stats->Add(lostContigs => 1);
    }
    $bin->add_prots(@$rolesL);
    if (! @$rolesL) {
        $stats->Add(uniFreeContigs => 1);
    }
    # Get the next line.
    $fields = $loader->GetLine(contig => $fh);
}
# Now all the bins are created, but we lack the coverage information. This is computed from the vector file.
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
# Open the output file.
open(my $oh, ">", $outputFile) || die "Could not open contig output file: $!";
# Loop through the contigs.
for my $contig (keys %contigs) {
    my $bin = $contigs{$contig};
    $bin->WriteContig($oh);
    $stats->Add(contigBinWritten => 1);
}
# Close the output file.
close $oh;
# All done.
print "All done.\n" . $stats->Show();
