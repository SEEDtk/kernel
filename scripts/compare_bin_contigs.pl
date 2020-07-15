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
use Stats;
use GEO;
use EvalCon;
use File::Copy::Recursive;

=head1 Compare Contig Distribution Between Two Binning Directories

    compare_bin_contigs.pl [ options ] sourceDir targetDir outDir

This script compares two directories of GTOs which were binned from the same sample.  For each bin in the source directory that has contigs in
common with the second directory, the number of base pairs and contigs in common will be output.

Two output files will be produced in the output directory.  C<summary.tbl> will list the name, quality, DNA size, and contig count of each bin.
C<detail.tbl> will list the DNA letters and contig count in common for each bin pair.  C<sankey.tbl> will contain the SankeyMatic
input for creating a diagram of the comparison on L<http://sankeymatic.com/build/>.

=head2 Parameters

The positional parameters are the names of the two binning directories and the output directory.  The source directory is first, and the target directory is second.

Progress messages will be written to the standard output.

The command-line options are as follows:

=over 4

=item recursive

Process all subdirectories of each input directory instead of the input directories themselves.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('sourceDir targetDir outDir',
        ['recursive|R', "process all subdirectories of the specified directories"]
        );
# Create a statistics object.
my $stats = Stats->new();
# Get the directories.
my ($sourceDir, $targetDir, $outDir) = @ARGV;
if (! $sourceDir) {
    die "No source directory specified.";
} elsif (! -d $sourceDir) {
    die "Source directory $sourceDir not found or invalid.";
} elsif (! $targetDir) {
    die "No target directory specified.";
} elsif (! -d $targetDir) {
    die "Target directory $targetDir not found or invalid.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir: $!";
} else {
    print "Output files will be written to $outDir.\n";
}
# Read the role files.
print "Loading role table.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
# Determine the directories to process.
my @samples;
if (! $opt->recursive) {
    push @samples, [$sourceDir, $targetDir, $outDir];
} else {
    print "Finding subdirectories to process.\n";
    opendir(my $dh, $sourceDir) || die "Could not open $sourceDir: $!";
    my @names = grep { (substr($_,0,1) ne '.') && -d "$sourceDir/$_" && -d "$targetDir/$_" } readdir $dh;
    print scalar(@names) . " subdirectories found.\n";
    push @samples, map { ["$sourceDir/$_", "$targetDir/$_", "$outDir/$_"] } @names;
}
print "Processing samples.\n";
for my $sample (@samples) {
    process_sample(@$sample);
}

print "All done.\n" . $stats->Show();

sub process_sample {
    my ($sourceDir, $targetDir, $outDir) = @_;
    if (! -d $outDir) {
        print "Creating output directory $outDir.\n";
        File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir: $!";
    }
    # Start the output files.
    print "Comparing $sourceDir to $targetDir in $outDir.\n";
    open(my $sh, '>', "$outDir/summary.tbl") || die "Could not open summary.tbl: $!";
    print $sh "sample_type\tbin_id\tbin_name\ttaxon_id\tgood_flag\tcompleteness\tcontamination\tconsistency\tdna_size\tcontig_count\n";
    open(my $oh, '>', "$outDir/detail.tbl") || die "Could not open detail.tbl: $!";
    print $oh "source_id\ttarget_id\tcontigs\tdna\n";
    open(my $kh, '>', "$outDir/sankey.tbl") || die "Could not open sankey.tbl: $!";
    # This hash maps each contig ID to its target bin.
    my %tContigs;
    # This hash tracks unbinned targets.
    my %tMissing;
    # This list contains all the target bin IDs.
    my @targets;
    # Get all the genomes in the target bin.
    my $geos = get_genomes($targetDir);
    # Loop through them, loading the contig hash.
    for my $geo (@$geos) {
        my $binID = $geo->id;
        push @targets, $binID;
        print "Scanning contigs in $binID.\n";
        my $contigH = $geo->contigHash;
        for my $contig (keys %$contigH) {
            $tContigs{$contig} = $binID;
            $stats->Add(targetContig => 1);
            $tMissing{$contig} = $geo->contigLen($contig);
        }
        # Output this genome's summary.
        write_genome($sh, "target", $geo);
        $stats->Add(targetBin => 1);
    }
    # Add a dummy target for no commonality.
    push @targets, "unbinned";
    # Now get the genomes in the source bin.
    $geos = get_genomes($sourceDir);
    # Loop through them, matching contigs.
    for my $geo (@$geos) {
        # This hash will track [dnaSize, contigs] for each target sample.
        my %targets = map { $_ => [0, 0] } @targets;
        my $contigH = $geo->contigHash;
        for my $contig (keys %$contigH) {
            $stats->Add(sourceContig => 1);
            my $len = $geo->contigLen($contig);
            my $target = $tContigs{$contig};
            if (! $target) {
                $stats->Add(sourceUnbinned => 1);
                $target = "unbinned";
            } else {
                $tMissing{$contig} = 0;
            }
            my $targetData = $targets{$target};
            $targetData->[0] += $len;
            $targetData->[1]++;
        }
        # Write the summary for this source bin.
        write_genome($sh, "source", $geo);
        $stats->Add(sourceBin => 1);
        # Write the details for this source bin.
        write_details($oh, $kh, $geo->id, \@targets, \%targets);
    }
    # Write the details for the unbinned targets.
    my %targets = map { $_ => [0, 0] } @targets;
    for my $contig (keys %tMissing) {
        my $len = $tMissing{$contig};
        if ($len > 0) {
            my $target = $tContigs{$contig};
            my $targetData = $targets{$target};
            $targetData->[0] += $len;
            $targetData->[1]++;
            $stats->Add(targetUnbinned => 1);
        }
    }
    write_details($oh, $kh, "unused", \@targets, \%targets);
    close $oh;
    close $sh;
    close $kh;
}


## Write detail lines.
sub write_details {
    my ($oh, $kh, $source, $targetL, $targetH) = @_;
    for my $target (@$targetL) {
        my $targetData = $targetH->{$target};
        my ($dnaSize, $contigCount) = @$targetData;
        if ($contigCount > 0) {
            print $oh join("\t", $source, $target, $contigCount, $dnaSize) . "\n";
            print $kh "$source [$contigCount] $target\n";
        }
    }
}
## Get all the GEOs from the specified directory.
sub get_genomes {
    my ($dir) = @_;
    print "Scanning $dir.\n";
    opendir(my $dh, $dir) || die "Could not open $dir: $!";
    my @gtoFiles = map { "$dir/$_" } grep { $_ =~ /\.gto$/ } readdir $dh;
    $stats->Add(bins => scalar @gtoFiles);
    my $gHash = GEO->CreateFromGtoFiles(\@gtoFiles, detail => 1, roleHashes => [$nMap, $cMap], logH => \*STDOUT, stats => $stats);
    my @genomes = sort keys %$gHash;
    return [ map { $gHash->{$_} } @genomes];
}

## Output summary data for a genome.
sub write_genome {
    my ($sh, $type, $geo) = @_;
    my $contigCount = $geo->contigCount;
    my $quality = $geo->quality;
    my $dnaLen = $quality->{metrics}{totlen};
    my $completeness = $quality->{complete};
    my $contamination = $quality->{contam};
    my $consistency = $quality->{fine_consis};
    my $goodFlag = ($geo->is_good ? 'Y' : '');
    my $name = $geo->name;
    my $id = $geo->id;
    my $taxon = $geo->taxon;
    print $sh join("\t", $type, $id, $name, $taxon, $goodFlag, $completeness, $contamination, $consistency, $dnaLen, $contigCount) . "\n";
}
