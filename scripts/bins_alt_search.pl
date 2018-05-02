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
use BasicLocation;
use Loader;
use Bin::Blast;

=head1 Search for Alternate Reference Genomes in Bin Jobs

    bins_alt_search.pl [ options ] binDir protDB

This script runs through the binning jobs in a specified directory and searches for alternative reference genomes in a specified
protein database. If one is found, the name of the bin job and the ID of the genome will be written to the standard output.
The genome can be tested for quality, and if it is good, added to the master seed protein database. The job can be rerun with
the new database to see if it improves the bins.

Status messages will be written to the standard error output.

=head2 Parameters

The positional parameters are the name of the bin job directory and the name of the BLAST database for the SEED proteins.
The command-line options are as follows.

=over 4

=item maxE

The maximum acceptable E-value. The default is C<1e-10>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=back

=cut

# Turn off buffering for stdout.
$| = 1;
# Create the loader object.
my $loader = Loader->new();
my $stats = $loader->stats;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir protDB',
                ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-10 }],
                ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
                ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
        );
# Check the parameters.
my ($binDir, $protDB) = @ARGV;
if (! $binDir) {
    die "No bin job directory specified.";
} elsif (! -d $binDir) {
    die "Bin job directory $binDir invalid or not found.";
} elsif (! $protDB) {
    die "No seed protein database specified.";
} elsif (! -s $protDB) {
    die "Seed protein database $protDB missing or empty.";
}
# Get the options.
my $maxE = $opt->maxe;
# Get the list of bin jobs.
print STDERR "Scanning $binDir.\n";
opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
my @binJobs = grep { -s "$binDir/$_/bins.found.tbl" } readdir $dh;
closedir $dh; undef $dh;
my $totJobs = scalar @binJobs;
my $countJobs = 0;
print STDERR "$totJobs bin jobs found in $binDir.\n";
# Loop through the bin jobs.
for my $jobDir (@binJobs) {
    my $workDir = "$binDir/$jobDir";
    $stats->Add(jobsProcessed => 1);
    $countJobs++;
    print STDERR "Processing $jobDir ($countJobs of $totJobs).\n";
    my $reducedFastaFile = "$workDir/reduced.fasta";
    my $binsFoundFile = "$workDir/bins.found.tbl";
    # Create the BLAST utility object.
    my $blaster = Bin::Blast->new(undef, $workDir, $reducedFastaFile,
        maxE => $maxE, minlen => $opt->minlen, gap => $opt->gap);
    # We will put the seed protein match locations in here.
    my %matches;
    print STDERR "Reading bin definitions from $binsFoundFile.\n";
    my $ih = $loader->OpenFile(binsFound => $binsFoundFile);
    while (my $binFields = $loader->GetLine(binsFound => $ih)) {
        my ($contig, $begin, $dir, $len) = @$binFields;
        $matches{$contig} = BasicLocation->new($contig, $begin, $dir, $len);
        $stats->Add(binsProcessed => 1);
    }
    close $ih; undef $ih;
    print STDERR "Reading old reference genomes.\n";
    my %oldGenomes;
    $ih = $loader->OpenFile(refGenomes => "$workDir/ref.genomes.scores.tbl");
    while (my $refData = $loader->GetLine(refGenomes => $ih)) {
        my (undef, $genome, undef, $name) = @$refData;
        $oldGenomes{$genome} = $name;
        $stats->Add(refsProcessed => 1);
    }
    close $ih; undef $ih;
    print STDERR "Searching for new reference genomes.\n";
    # Create a hash mapping each contig ID to the DNA sequence representing the hit. We do this by reading
    # the sample FASTA file and applying the matches hash.
    my $seqHash = $loader->GetDNA(\%matches, $reducedFastaFile);
    # Now BLAST against the seed protein database.
    my $contigHash = $blaster->MatchProteins($seqHash, 'PhenTrnaSyntAlph', 1, $maxE, db => $protDB, type => 'dna');
    # Loop through the new hits.
    my $new = 0;
    for my $contig (sort keys %$contigHash) {
        for my $genomeData (@{$contigHash->{$contig}}) {
            my ($genome, $score, $name) = @$genomeData;
            if (! $oldGenomes{$genome}) {
                print join("\t", $jobDir, $genome, $name) . "\n";
                $new++;
                $stats->Add(newGenomeFound => 1);
            }
        }
    }
    if ($new) {
        print STDERR "$new new reference genome(s) found.\n";
    }
}
print STDERR "All done.\n" . $stats->Show();