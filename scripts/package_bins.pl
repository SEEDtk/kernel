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
use Bin;
use GenomeTypeObject;
use File::Copy::Recursive;
use Stats;

=head1 Create Genome Packages

    package_bins.pl [ options ] sampleDir packageDir

This script will create genome packages from the bins in a single sample directory. Each genome package is a directory with the
same name as the bin's genome ID and contains the following files

=over 4

=item bin.gto

A L<GenomeTypeObject> for the genome formed from the bin.

=item bin.fa

A FASTA file containing the bin's contigs.

=item data.tbl

A tab-delimited file of key/value pairs, with the following keys.

=over 8

=item Genome Name

Name of the bin genome

=item Sample Name

Name of the sample from which the bin came.

=item Bin Number

The ID number of the bin in the sample.

=item Contigs

The number of contigs in the bin.

=item Base pairs

The number of DNA base pairs in the bin's contigs.

=item Ref Genome

The ID of the closest reference genome.

=item Ref Name

The name of the closest reference genome.

=back

=back

=head2 Parameters

The two positional parameters are the full name of the sample directory and the full name of the directory to contain the output
genome packages. If the genome package already exists, it will not be overwritten.

The command-line options are the following.

=over 4

=item force

If specified, existing genome packages will be overwritten.

=item recursive

If specified, then the sample directory is treated as a directory of sample directories, and all samples therein will be processed.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sampleDir packageDir',
        ['force', "overwrite existing packages"],
        ['recursive', "process a directory of samples"],
        );
my $stats = Stats->new();
# Get the directories.
my ($sampleDir, $packageDir) = @ARGV;
if (! $sampleDir) {
    die "No sample directory specified.";
} elsif (! -d $sampleDir) {
    die "Invalid sample directory $sampleDir.";
}
if (! $packageDir) {
    die "No package directory specified.";
} elsif (! -d $packageDir) {
    die "Invalid package directory $packageDir";
}
my @samples;
if (! $opt->recursive) {
    push @samples, $sampleDir;
} else {
    opendir(my $dh, $sampleDir) || die "Could not open sample directory: $!";
    push @samples, map { "$sampleDir/$_" } grep { substr($_, 0, 1) ne '.' && -d "$sampleDir/$_" } readdir $dh;
}
for my $sample (@samples) {
    $stats->Add(samples => 1);
    if (! -s "$sample/bins.rast.json") {
        print "Skipping $sample: incomplete.\n";
        $stats->Add(sampleSkipped => 1);
    } elsif (! -s "$sample/bin1.gto") {
        print "Skipping $sample: no bins.\n";
        $stats->Add(sampleEmpty => 1);
    } else {
        # Get the sample name.
        my @pieces = split /[\\\/:]/, $sample;
        my $sampleName = pop @pieces;
        print "Sample name is $sampleName.\n";
        # Read in the scores for the reference genomes. We use these to compute the closest genome to each bin.
        # The hash will map each genome ID to a [score, name] tuple.
        my %refGenomes;
        print "Processing reference genomes.\n";
        open(my $ih, '<', "$sample/ref.genomes.scores.tbl") || die "Could not open ref genomes score table: $!";
        while (! eof $ih) {
            my $line = <$ih>;
            if ($line =~ /(\d+\.\d+)\t([\d\.]+)\t(.+)$/) {
                $refGenomes{$1} = [$2, $3];
            }
        }
        # Now read in the bins.
        print "Reading bin file.\n";
        my $binList = Bin::ReadBins("$sample/bins.rast.json");
        # Create a map from contig IDs to bins.
        my %binMap;
        for my $bin (@$binList) {
            $binMap{$bin->contig1} = $bin;
        }
        # Now we process the bins sequentially. For each GTO/FA pair, we search the contigs to find the matching bin object.
        # The bin object produces the majority of the data.tbl values.
        my $binN = 1;
        while (-f "$sample/bin$binN.gto") {
            my $binName = "bin$binN";
            print "Processing $binName.\n";
            $stats->Add(bins => 1);
            my $gto = GenomeTypeObject->create_from_file("$sample/$binName.gto");
            # Get the ID and name of the genome.
            my $genomeID = $gto->{id};
            my $name = $gto->{scientific_name};
            # Compute the output directory name.
            my $genomeDir = "$packageDir/$genomeID";
            if (-d $genomeDir && ! $opt->force) {
                print "Package already exists for $binName: $genomeID.\n";
                $stats->Add(binsAlreadyFound => 1);
            } else {
                # Here the bin is new or we are forcing.
                if (! -d $genomeDir) {
                    print "Creating $genomeDir.\n";
                    mkdir $genomeDir;
                    $stats->Add(binProcessed => 1);
                } else {
                    print "Replacing $genomeDir.\n";
                    File::Copy::Recursive::pathempty($genomeDir);
                    $stats->Add(binReplaced => 1);
                }
                # Find the bin object. One of the contigs will identify the bin. We stop when we hit it.
                my $bin;
                my $contigs = $gto->{contigs};
                for my $contig (@$contigs) { last if $bin;
                    $bin = $binMap{$contig->{id}};
                }
                die "Bin object not found for $binName." if ! $bin;
                # Copy the main files.
                File::Copy::Recursive::fcopy("$sample/$binName.gto", "$genomeDir/bin.gto") ||
                    die "Error copying $binName.gto: $!";
                File::Copy::Recursive::fcopy("$sample/$binName.fa", "$genomeDir/bin.fa") ||
                    die "Error copying $binName.fa: $!";
                # We will compute the data table values in here.
                my %data;
                $data{'Genome Name'} = $name;
                $data{'Sample Name'} = $sampleName;
                $data{'Bin Number'} = $binN;
                $data{'Contigs'} = $bin->contigCount;
                $data{'Base pairs'} = $bin->len;
                # Find the closest reference genome.
                my @refs = $bin->refGenomes;
                my ($closest, $best) = ('', 0);
                for my $ref (@refs) {
                    my $score = $refGenomes{$ref}[0];
                    if ($score >= $best) {
                        $best = $score;
                        $closest = $ref;
                    }
                }
                $data{'Ref Genome'} = $closest;
                $data{'Ref Name'} = $refGenomes{$closest}[1];
                # Write the data file.
                open(my $oh, '>', "$genomeDir/data.tbl") || die "Could not open $binName data table file: $!";
                for my $key (sort keys %data) {
                    print $oh "$key\t$data{$key}\n";
                }
            }
            # Move to the next bin.
            $binN++;
        }
    }
}
print "All done.\n" . $stats->Show();