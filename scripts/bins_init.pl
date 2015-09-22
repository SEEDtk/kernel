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
use KmerRepGenomes;
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

=item lenFiter

Minimum contig length for a contig to be considered meaningful. Only meaningful contigs are retained. The default
is C<10000>.

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

=item workDir

Working directory to contain the intermediate files. If omitted, the sample directory is assumed.

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

=item maxrefs

The maximum number of reference genomes to keep for BLASTing purposes.

=back

=head2 Working Files

This script produces intermediate working files that are re-used if they already exist. Many files are output in the
working directory by L<Bin::Blast> and L<KmerRepGenomes>. These are independent of the contigs input, so the same working
directory can be used for multiple communities.

The C<ref.counts.tbl> file contains the IDs of the reference genomes to be used for the analysis of the sample.
This is stored in the sample directory.

The C<sample.fasta> file contains the sequences for the meaningful contigs in the sample.

The C<rejected.fasta> file contains the sequences for the non-meaningful contigs in the sample.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir outputFile',
                Shrub::script_options(),
                ['minsim=f',     'minimum percent identity score for BLAST results', { 'default' => 0.9 }],
                ['minlen=i',     'minimum length for an acceptable BLAST match', { 'default' => 150 }],
                ['tetra=s',      'tetranucleotide counting algorithm', { 'default' => 'dual' }],
                ['workDir|D=s',  'working directory for intermediate files'],
                ['refs=s',       'reference genome list file'],
                ['unifile=s',    'universal role file', { default => "$FIG_Config::global/uni_roles.tbl" }],
                ['reps=s',       'representative role file', { default => "$FIG_Config::global/rep_roles.tbl" }],
                ['force',        'force a rebuild of the BLAST and KMER databases'],
                ['minhits=i',    'minimum number of kmer hits for a genome to be considered useful', { default => 400 }],
                ['maxfound=i',   'maximum number of kmer occurrences before it is considered common', { default => 10 }],
                ['kmersize|k=i', 'kmer size (protein)', { default => 10 }],
                ['lenFilter=i',  'minimum contig length', { default => 10000 }],
                ['covgFilter=f',  'minimum contig mean coverage', { default => 10}],
                ['maxrefs=i',    'maximum number of reference genomes to use', { default => 10 }],
        );
# Turn off buffering for stdout.
$| = 1;
# Create the loader object.
my $loader = Loader->new();
# Get the statistics object.
my $stats = $loader->stats;
# Save the force-flag and the filter options.
my $force = $opt->force // 0;
my $lenFilter = $opt->lenfilter;
my $covgFilter = $opt->covgfilter;
# Check the file names.
my ($sampleDir, $outputFile) = @ARGV;
if (! $sampleDir) {
    die "Sample directory name missing.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
my ($contigFile, $vectorFile, $reducedFastaFile, $rejectedFastaFile) =
        map { "$sampleDir/$_" } qw(contigs.fasta output.contigs2reads.txt sample.fasta rejected.fasta);
if (! -f $contigFile) {
    die "Contig file $contigFile not found.";
} elsif (! -f $vectorFile) {
    die "Vector file $vectorFile not found.";
} elsif (! $outputFile) {
    die "Output file name missing.";
}
# Check the working directory.
my $workDir = $opt->workdir;
if (! $workDir) {
    $workDir = $sampleDir;
} elsif (! -d $workDir) {
    die "Invalid working directory $workDir.";
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
$fh = $loader->OpenFasta(contig => $contigFile);
$fields = $loader->GetLine(contig => $fh);
while (defined $fields) {
    my ($contig, undef, $seq) = @$fields;
    if (exists $contigs{$contig}) {
        print $ofh ">$contig\n$seq\n";
    } else {
        print $xfh ">$contig\n$seq\n";
    }
    $fields = $loader->GetLine(contig => $fh);
}
close $ofh;
close $xfh;
my $refsFile = "$sampleDir/ref.counts.tbl";
# Get the list of representative role IDs.
my $repRoles = $loader->GetNamesFromFile(repRoles => $opt->reps);
my $repCount = scalar @$repRoles;
print "$repCount representative roles found.\n";
# Create the KMER database.
print "Creating KMER database.\n";
my $kmers = KmerRepGenomes->new($shrub, "$workDir/rep_kmers.json", $repRoles,
        force => $force, kmerSize => $opt->kmersize, minHits => $opt->minhits,
        maxFound => $opt->maxfound);
# Compute the list of genomes to BLAST. First we find the relevant genomes.
print "Computing reference genomes.\n";
my $refGenomeH = $kmers->FindGenomes($reducedFastaFile);
my $closeCount = scalar keys %$refGenomeH;
$stats->Add(closeRefGenome => $closeCount);
print "$closeCount close genomes found.\n";
# Filter the complete list of genomes by the known representatives.
my @filtered;
for my $refGenome (@$refs) {
    my $count = $refGenomeH->{$refGenome};
    if ($count) {
        push @filtered, $refGenome;
        $stats->Add(keptRefGenome => 1);
    }
}
my @sortedRefs = sort { $refGenomeH->{$b}[0] <=> $refGenomeH->{$a}[0] } @filtered;
my @realRefs;
if (scalar(@sortedRefs) > $opt->maxrefs) {
    @realRefs = @sortedRefs[0 .. $opt->maxrefs];
} else {
    @realRefs = @sortedRefs;
}
# Prepare to write what we decided to keep.
open(my $orh, ">", $refsFile) || die "Could not open reference genome save file: $!";
for my $refGenome (@realRefs) {
    print $orh join("\t", $refGenome, @{$refGenomeH->{$refGenome}}) . "\n";
}
print scalar(@realRefs) . " reference genomes qualified.\n";
# Finally, we do the BLASTing.
my $blastStats = Bin::Blast::Process($shrub, $sampleDir, \@realRefs, \%contigs, $reducedFastaFile, $uniRoles,
        minsim => $opt->minsim, minlen => $opt->minlen);
$stats->Accumulate($blastStats);
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
