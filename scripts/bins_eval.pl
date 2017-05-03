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
use Shrub;
use Bin;
use Bin::Blast;
use Bin::Analyze;
use File::Copy::Recursive;
use Loader;

=head1 Evaluate Community Bins

    bins_eval.pl [ options ] binDir

This program computes a more realistic appraisal of the universal roles occurring in a bin. For each major bin, we take the
first reference genome and extract the proteins for its universal roles. These are BLASTed against the DNA for the bin.

=head2 Parameters

There is one positional parameter-- the directory containing the bin information. The json file containing the bin data structures
must be named C<bins.json> and the FASTA data for all the sample contigs must be in C<sample.fasta>. These are output by
L<bins_create.pl>.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item unifile

The name of a file containing the IDs of the universal roles. If omitted, the file C<uni_roles.tbl> in the C<Global>
subdirectory of the data directory will be used. The file is tab-delimited, with the role IDs in the first column,
occurrence counts in the second, and role names in the third.

=item maxE

The maximum acceptable E-value. The default is C<1e-30>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDirectory', Shrub::script_options(),
                ['unifile=s',      'universal role file', { default => "$FIG_Config::global/uni_roles.tbl" }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-50 }],
                ['minlen|p=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }]
        );
# Create the loader object and get the statistics.
my $loader = Loader->new();
my $stats = $loader->stats;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Read the bins.
my $binList;
my ($binDir) = @ARGV;
if (! $binDir) {
    die "A bins directory is required.";
} elsif (! -d $binDir) {
    die "Bin directory $binDir not found or invalid.";
}
# Compute the bins file name.
my $jsonFileName = "$binDir/bins.json";
print "Reading bins.\n";
$binList = Bin::ReadBins($jsonFileName);
# Create a directory hash that matches each contig to a bin. Only bins with a reference genome are
# kept.
my %contig2Bin;
my @goodBins;
for my $bin (@$binList) {
    $stats->Add(binFound => 1);
    my ($refGenome) = $bin->refGenomes;
    if ($refGenome) {
        $stats->Add(binkept => 1);
        for my $contig ($bin->contigs) {
            $contig2Bin{$contig} = scalar(@goodBins);
            $stats->Add(contigInBin => 1);
        }
        push @goodBins, $bin;
    }
}
# Create a FASTA output file for each bin.
my @binHandles;
for (my $id = 0; $id < scalar(@goodBins); $id++) {
    open(my $bh, ">$binDir/bin$id.fa") || die "Could not create FASTA file for bin $id";
    push @binHandles, $bh;
}
# Get a BLAST object.
my $blaster = Bin::Blast->new($shrub, $binDir, "$binDir/sample.fasta", uniRoles => $opt->unifile,
        maxE => $opt->maxe, minlen => $opt->minlen, gap => $opt->gap);
# Now we read the contigs from the sample.
print "Sorting contigs.\n";
my $ih = $loader->OpenFasta(sampleContigs => "$binDir/sample.fasta");
while (my $fields = $loader->GetLine(sampleContigs => $ih)) {
    my ($contig, $comment, $dna) = @$fields;
    my $contigBin = $contig2Bin{$contig};
    if (defined $contigBin) {
        # Here this contig belongs in a bin. Write it to the bin's FASTA file.
        my $bh = $binHandles[$contigBin];
        print $bh ">$contig $comment\n$dna\n";
        $stats->Add(contigBinOut => 1);
    }
}
# Close the bin FASTA files.
for my $bh (@binHandles) {
    close $bh;
}
# Loop through the bins we are interested in.
for (my $id = 0; $id < scalar(@goodBins); $id++) {
    my $binFasta = "$binDir/bin$id.fa";
    my $bin = $goodBins[$id];
    print "Processing bin $id.\n";
    # Get the reference genome. We know one exists, because that is the definition of a good bin.
    my ($genome) = $bin->refGenomes;
    # Find the universal roles that match bin DNA.
    print "Searching for best unirole hits in $binFasta.\n";
    my $protHitHash = $blaster->GoodProteinHits($genome, $binFasta);
    my $protsFound = scalar(keys %$protHitHash);
    print "$protsFound uniroles found.\n";
    $stats->Add(proteinsFound => $protsFound);
    # Clear the bin's universal protein list.
    $bin->replace_prots();
    # Store the new counts.
    for my $prot (keys %$protHitHash) {
        my $hitList = $protHitHash->{$prot};
        $bin->incr_prot($prot => scalar @$hitList);
    }
}
print "Writing new bins.\n";
open(my $oh, ">$binDir/bins.new.json") || die "Could not open bins output file: $!";
for my $bin (@$binList) {
    $bin->Write($oh);
}
# Get the universal role information.
print "Producing report.\n";
my $uniRoles = $blaster->uni_hash;
my $totUnis = $blaster->uni_total;
# Produce the bin report.
open(my $rh, ">$binDir/bins.report.good.txt") || die "Could not open report file: $!";
my $analyzer = Bin::Analyze->new(totUnis => $totUnis, minUnis => (0.8 * $totUnis));
$analyzer->BinReport($rh, $uniRoles, \@goodBins);
close $rh;
print "All done.\n" . $stats->Show();