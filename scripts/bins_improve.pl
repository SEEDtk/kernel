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
use RASTlib;
use FastA;
use GEO;
use EvalCon;
use P3DataAPI;
use Stats;
use SeedUtils;

=head1 Attempt to Improve Bins by Removing Contamination

    bins_improve.pl [ options ] binDir

This script will search the bins in a binning directory.  Bins with a high completeness but also a high contamination will be examined
for contigs that can be removed, and then will be re-submitted for annotation.

A deep-evaluation L<GEO> with reference genomes will be used to identify which contigs have no good features.  A feature is good if it
implements a good role (one whose predictions match the actual annotations) or if it is close to a feature with the same role in the reference
genome.

The C<bins.report.txt> file will be read to get a list of the bins and their reference genomes.  Then the individual C<bin>I<N>C<.gto> files
will be read and analyzed for quality information.  If the bin's completeness and contamination are both high, the genomes and reference
genomes will be read in as GEO objects and the bad contigs identified.  A new C<bin>I<N>C<A.fa> file will be produced and
submitted for annotation using L<RASTlib>.

=head2 Parameters

The positional parameter is the name of the binning directory.

The command-line options are as follows.

=over 4

=item recursive

If specified, then the input directory is presumed to be a master directory whose subdirectories should all be processed.

=item noindex

If specified, the new genomes will not be indexed in PATRIC.

=item force

If specified, bins that have already been improved will be re-processed.  Normally, if a C<bin>I<N>C<A.gto> file exists, the bin
will be skipped.

=item sleep

Sleep interval for RAST status polling. The default is C<60>.

=item minComplete

The minimum completeness score for a bin to be considered for improvement.  The default is C<90>.

=item minContam

The minimum contamination score for a bin to be considered for improvement.  The default is C<10>.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('binDir',
            ['recursive', 'process subdirectories'],
            ['noindex', 'do not integrate new genomes into the PATRIC index'],
            ['force', 're-process bins that have already been improved'],
            ["sleep=i", "sleep interval for RAST status polling", { default => 60 }],
            ["minComplete|complete|M=i", 'minimum completeness score', { default => 90 }],
            ["minContam|contam|m=i", 'minimum contamination score', { default => 10 }]
        );
# Create the support objects.
my $stats = Stats->new();
my $p3 = P3DataAPI->new();
# Get the options.
my $noIndex = $opt->noindex // 0;
my $force = $opt->force;
my $sleep = $opt->sleep;
my $minComplete = $opt->mincomplete;
my $minContam = $opt->mincontam;
# Compute the directories to process.
my @binDirs;
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "$binDir is missing or invalid.";
} elsif (! $opt->recursive) {
    push @binDirs, $binDir;
    print "Input directory is $binDir.\n";
} else {
    print "Scanning $binDir.\n";
    opendir(my $dh, $binDir) || die "Could not open directory $binDir: $!";
    @binDirs = map { "$binDir/$_" } grep { -s "$binDir/$_/bins.report.txt" } readdir $dh;
    closedir $dh;
    print scalar(@binDirs) . " input directories found in $binDir.\n";
}
# Read the role files.
print "Reading role files.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::global/roles.in.subsystems", $stats);
# Create the GEO options.
my %gOptions = (roleHashes => [$nMap, $cMap], detail => 2, p3 => $p3, stats => $stats);
# Create the RAST options.
my %rOptions = (noIndex => $noIndex, 'sleep' => $sleep);
# Loop through the binning directories.
for my $sample (@binDirs) {
    print "Processing directory $sample.\n";
    open(my $ih, '<', "$sample/bins.report.txt") || die "Could not open bins.report.txt: $!";
    # We scan this file to pull out the reference genomes for each bin.
    my %refs;
    my $bin = 0;
    while (! eof $ih) {
        my $line = <$ih>;
        $stats->Add(binReportLine => 1);
        if ($line =~ /^BIN\s+(\d+)/) {
            # Here we have a bin header.
            $bin = $1;
            $stats->Add(binHeader => 1);
        } elsif ($line =~ /^\s+(\d+\.\d+):\s+(.+)/) {
            print "Bin $bin has reference $1: $2.\n";
            push @{$refs{$bin}}, $1;
            $stats->Add(binRefLine => 1);
        }
    }
    # Now we have a list of bins to process.  We read the GTO of each bin and extract the quality data.
    for my $bin (sort { $a <=> $b } keys %refs) {
        my $gtoFile = "$sample/bin$bin.gto";
        my $fastaFile = "$sample/bin$bin.fa";
        my $gtoNew = "$sample/bin${bin}A.gto";
        my $fastaNew = "$sample/bin${bin}A.fa";
        if (-s $gtoNew && ! $force) {
            print "$sample bin $bin already processed-- skipping.\n";
            $stats->Add(binAlreadyProcessed => 1);
        } else {
            # Create the GEO.
            print "Creating GEO for $sample bin $bin.\n";
            my $geo = GEO->CreateFromGto($gtoFile, %gOptions);
            # Get the quality data.
            my ($coarse, $fine, $complete, $contam, $group) = $geo->scores;
            unless ($complete >= $minComplete && $contam >= $minContam) {
                $stats->Add(binNotQualified => 1);
                print "$sample bin $bin has completeness $complete and contamination $contam-- skipping.\n";
            } else {
                print "Loading reference genomes for $sample bin $bin.\n";
                my $gHash = GEO->CreateFromPatric($refs{$bin}, %gOptions);
                for my $ref (keys %$gHash) {
                    $geo->AddRefGenome($gHash->{$ref});
                }
                # Do a full quality analysis.
                $geo->AnalyzeQualityData();
                # Get the contig list.  We will put the bad contig IDs in here.
                print "Analyzing contigs for $sample bin $bin.\n";
                my %bads;
                my $contigs = $geo->{quality}{contigs};
                for my $contig (keys %$contigs) {
                    $stats->Add(contigChecked => 1);
                    if (! $contigs->{$contig}[0]) {
                        $stats->Add(badContig => 1);
                        $bads{$contig} = 1;
                    }
                }
                # Now we set up for the RASTlib call.  This involves computing the parameters and reading the contigs.
                # We will also create a new FASTA file to contain the contigs.
                my $taxon = $geo->taxon;
                my $name = $geo->name . " cleaned";
                my $gc = $geo->gc;
                my $domain = substr($geo->domain, 0, 1);
                # Now read the contigs and create a new FASTA file.
                my @triples;
                print "Reading contigs for $sample bin $bin.\n";
                open(my $oh, '>', $fastaNew) || die "Could not open FASTA output file: $!";
                my $fh = FastA->new($fastaFile);
                while ($fh->next) {
                    $stats->Add(fastaIn => 1);
                    my $id = $fh->id;
                    my $seq = $fh->left;
                    my $len = length $seq;
                    $stats->Add(dnaIn => $len);
                    if ($bads{$id}) {
                        $stats->Add(dnaSkipped => $len);
                    } else {
                        $fh->Write($oh);
                        push @triples, [$id, '', $seq];
                        $stats->Add(dnaOut => $len);
                        $stats->Add(fastaOut => 1);
                    }
                }
                undef $fh;
                close $oh;
                # Submit this bin to RAST.
                print "Submitting $fastaNew to RAST.\n";
                my $gto = RASTlib::Annotate(\@triples, $taxon, $name, %rOptions, geneticCode => $gc, domain => $domain);
                # Look at the results.  We use detail level 0 because we don't need details any more.
                $geo = GEO->CreateFromGto($gto, %gOptions, detail => 0);
                my $id2 = $geo->id;
                my $name2 = $geo->name;
                my $status = 'still bad';
                if ($geo->is_good) {
                    $status = "good";
                    $stats->Add(binImproved => 1);
                }
                my ($coarse2, $fine2, $complete2, $contam2, $group2) = $geo->scores;
                print "New genome is $id2: $name2 and is $status.\n";
                print "NEW SCORES: complete $complete to $complete2, fine $fine to $fine2, contamination $contam to $contam2.\n";
                SeedUtils::write_encoded_object($gto, $gtoNew);
                print "$gtoNew saved.\n";
            }
        }
    }
}
print "All done.\n" . $stats->Show();