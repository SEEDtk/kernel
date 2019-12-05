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

=head1 Close-Genome Annotation Stress Test

    Stress_Close_Anno.pl [options] testDir

This script stress-tests the close-genome annotation process.  The test directory should contain a C<.tbl> file for
each target genome to be tested.  The file will be tab-delimited with headers, with genome IDs in the first column,
sorted from closest to furthest.  For each target genome, we want to test the closest genomes, the furthest genomes,
the two closest, the one furthest, and all the genomes.  For each test, we create a working directory and a control
file and produce the output in there.  The output will include a GTO, a result table, a result log, a comparison table,
and a comparison log.

To pre-process a genome, we download a GTO for the genome and remove all the features, then download protein FASTA
files for all the genomes in its close-genome table.  Then for each test, we create the control file and the
working directory, run L<Close_Anno.pl>, and finally run L<p3x-genome-function-file.pl>.

The control files created will have extra fields used to analyze the results.  In addition to the FASTA file name and
the genetic code, there will be a column for the GC content and the normalized similiarity with the target genome.
These are used to create a summary of the results.

=head2 Parameters

The positional parameter is the name of the testing directory.

The command-line options are as follows:

=over 4

=item build

If specified, then the input files should be built and the genomes pre-processed.  Otherwise, the test directories are presumed to already exist
and the tests are simply run from inside them.  The parameter should be the name of a file containing the IDs of the target genomes in the first
column.  The incoming genomes must be representative genomes found in the file C<rep100.sorted.tbl> in the SEEDtk global directory.  This
file contains all the good PATRIC genomes sorted by similarity to its representative.

=item missing

Only run a test if it has not already been run.

=item compare

Only run the comparisons.

=back

=cut

use strict;
use FIG_Config;
use P3DataAPI;
use P3Utils;
use GenomeTypeObject;
use IPC::Run3;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('testDir',
        ['build=s', 'build the test directories and files from the specified input genomes'],
        ['compare', 'rerun the comparisons only'],
        ['missing', 'only run tests that have not been run yet']);
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the testing directory.
my ($testDir) = @ARGV;
if (! $testDir) {
    die "No testing directory specified.";
} elsif (! -d $testDir) {
    die "Invalid or missing testing directory $testDir.";
}
if ($opt->build) {
    open(my $ih, '<', $opt->build) || die "Could not open build file: $!";
    BuildTestDirectory($ih, $testDir);
}
# Get the test genome directories.
opendir(my $dh, $testDir) || die "Could not open $testDir: $!";
my @gDirs = sort grep { $_ =~ /^\d+\.\d+\./ && -f "$testDir/$_/control.tbl" } readdir $dh;
print scalar(@gDirs) . " test directories found.\n";
# Make the test directory current.
chdir "$testDir";
# Loop through the test genome directories.
for my $gDir (@gDirs) {
    RunTest($gDir);
}
# Now build the analysis file.
print "Analyzing results.\n";
my @stats = qw(test_name patricPeg newPeg samePeg sameLength sameFunction closeLength extraOrf lostOrf patricLonger newLonger duration pct_found pct_close pct_same similarity neighbors gc_content);
open(my $oh, '>', "analysis.tbl") || die "Could not open analysis file: $!";
P3Utils::print_cols(\@stats, oh => $oh);
for my $gDir (@gDirs) {
    print "Analyzing $gDir.\n";
    my %columns = map {$_ => 0} @stats;
    $columns{test_name} = $gDir;
    open(my $ih, '<', "$gDir/DONE") || die "Could not open completion file for $gDir: $!";
    $columns{duration} = <$ih>;
    close $ih; undef $ih;
    open($ih, '<', "$gDir/control.tbl") || die "Could not open control file for $gDir: $!";
    my ($score, $content, $count) = (0, 0, 0);
    while (! eof $ih) {
        my ($genomeID, $gc, $scoreX, $contentX) = P3Utils::get_fields($ih);
        $score += $scoreX;
        $content += $contentX;
        $count++;
    }
    $columns{gc_content} = $content / $count;
    $columns{similarity} = $score / $count;
    $columns{neighbors} = $count;
    close $ih; undef $ih;
    open($ih, '<', "$gDir/compare.log") || die "Could not open comparison file for $gDir: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\S+)\s+(\d+)$/) {
            if (exists $columns{$1}) {
                $columns{$1} = $2;
            }
        }
    }
    close $ih;
    my $orfs = $columns{patricPeg};
    $columns{pct_found} = ($orfs - $columns{lostOrf}) * 100 / $orfs;
    $columns{pct_close} = $columns{closeLength} * 100 / $orfs;
    $columns{pct_same} = $columns{samePeg} * 100 / $orfs;
    P3Utils::print_cols([map { $columns{$_} } @stats], oh => $oh);
}
close $oh;
print "All done.\n";

## Create the test directory.  For each genome, we need a contigs-only GTO file in the test directory, the protein FASTA files for its
## neighbors, and 3 or more subdirectories set up to run the tests.
sub BuildTestDirectory {
    my ($ih, $testDir) = @_;
    if (-d $testDir) {
        print "Clearing test directory.\n";
        File::Copy::Recursive::pathempty($testDir) || die "Could not erase $testDir: $!";
    } else {
        print "Creating test directory.\n";
        File::Copy::Recursive::pathmk($testDir) || die "Could not create $testDir: $!";
    }
    # Now read in the target genome IDs.
    print "Reading target genomes.\n";
    my %genomes;
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /(\d+.\d+)/) {
            $genomes{$1} = [];
        }
    }
    my $gCount = scalar keys %genomes;
    print "$gCount genome IDs found.\n";
    # Now we read the rep-genomes file to find the neighbors.  For each genome, there will be a sub-list of [neighbor, score]
    # pairs.
    my %normalizer;
    open(my $rh, '<', "$FIG_Config::global/rep100.sorted.tbl") || die "Could not open rep100.sorted.tbl: $!";
    # Skip the header.
    my $line = <$rh>;
    while (! eof $rh) {
        my ($genome, $name, $rep, $score) = P3Utils::get_fields($rh);
        if ($genomes{$genome}) {
            print "$genome ($name) found.\n";
            # Save the normalization factor for this genome.
            $normalizer{$genome} = $score;
        } elsif ($genomes{$rep}) {
            # Here we have a possible neighbor.
            my $subList = $genomes{$rep};
            if (scalar(@$subList) < 10) {
                push @$subList, [$genome, $score];
                print "$genome ($name) is a neighbor of $rep.\n";
            }
        }
    }
    close $rh;
    # We have the neighbor lists.  Now we process each genome and its neighbors.
    for my $genome (sort keys %genomes) {
        CreateSourceGto($genome);
        # Get the neighbor list for this genome.
        my $nList = $genomes{$genome};
        my $nDenom = $normalizer{$genome};
        # Create the neighbor FASTA files.  First we need the genetic code and GC content.
        my @nCouplets = map { [$_->[0], [$_->[0], $_->[1]/$nDenom]] } @$nList;
        print "Retrieving $genome neighbor data.\n";
        my $nData = P3Utils::get_data_batch($p3, genome => [], ['gc_content', 'genetic_code'], \@nCouplets, 'genome_id');
        my $nCount = scalar @$nData;
        print "$nCount neighbors found.\n";
        die "Two few neighbors for $genome.\n" if $nCount < 4;
        # Now we create the files themselves.
        for my $nTuple (@$nData) {
            my $nID = $nTuple->[0];
            print "Downloading $nID proteins.\n";
            my $proteins = P3Utils::get_data($p3, feature => [['eq', 'genome_id', $nID], ['eq', 'feature_type', 'CDS']], ['patric_id', 'product', 'aa_sequence']);
            open(my $fh, '>', "$testDir/$nID.faa") || die "Could not open $nID.faa: $!";
            my ($kept, $skipped) = (0, 0);
            for my $protein (@$proteins) {
                my $aaSeq = uc $protein->[2];
                if ($protein->[1] && $aaSeq && $aaSeq =~ /^[GALMFWKQESPVICYHRNDT]+$/) {
                    print $fh ">$protein->[0] $protein->[1]\n$aaSeq\n";
                    $kept++;
                } else {
                    print "Invalid protein $protein->[0] skipped.\n";
                    $skipped++;
                }
            }
            print "$kept proteins kept, $skipped skipped.\n";
        }
        # Finally, we create the testing sub-directories.
        CreateTestDir($genome, $nData, 0, closest => $nCount);
        if ($nCount > 5) {
            CreateTestDir($genome, $nData, 0, closest => 5);
            CreateTestDir($genome, $nData, $nCount - 5, furthest => 5);
        }
        CreateTestDir($genome, $nData, 0, closest => 2);
        CreateTestDir($genome, $nData, $nCount - 1, furthest => 1);
    }
}

# Create a single testing subdirectory.
sub CreateTestDir {
    my ($genome, $nData, $n0, $type, $nCount) = @_;
    # Create the subdirectory.
    my $testDirName = "$testDir/$genome.$type$nCount";
    print "Creating $testDirName.\n";
    File::Copy::Recursive::pathmk($testDirName) || die "Could not create $testDirName: $!";
    open(my $oh, '>', "$testDirName/control.tbl") || die "Could not open control file: $!";
    my $n2 = $n0 + $nCount;
    for (my $i = $n0; $i < $n2; $i++) {
        my ($genome, $score, $gc_content, $gc) = @{$nData->[$i]};
        print $oh "$genome.faa\t$gc\t$score\t$gc_content\n";
    }
    close $oh;
}

# Create the target genome's GTO in the main test directory.
sub CreateSourceGto {
    my ($genomeID) = @_;
    my $gtoName = "$testDir/$genomeID.gto";
    if (-s $gtoName) {
        print "$genomeID already downloaded.\n";
    } else {
        print "Downloading $genomeID.\n";
        my $gto = $p3->gto_of($genomeID);
        die "$genomeID not found." if (! $gto);
        print "Preparing GTO.\n";
        $gto->{features} = [];
        delete $gto->{quality} if $gto->{quality};
        $gto->{subsystems} = [];
        $gto->destroy_to_file($gtoName);
        print "$genomeID is $gto->{scientific_name} with GC content $gto->{gc_content} from $gto->{domain}.\n";
    }
}

## Run a test in a single test directory.
sub RunTest {
    my ($testDirName) = @_;
    # Get the genome ID.
    my ($genomeID) = ($testDirName =~ /^(\d+\.\d+)/);
    # Compute the commands.
    my $closeAnno = "Close_Anno --log $testDirName/results.log --gto $testDirName/results.gto $genomeID.gto $testDirName/control.tbl $testDirName";
    my $analyze = "p3x-genome-function-file --log $testDirName/compare.log --input $testDirName/results.tbl $genomeID";
    # Check to see if we need to run this test.
    my $done = -f "$testDirName/DONE";
    if ($opt->missing && $done) {
        print "Skipping $testDirName for $genomeID:  already found.\n";
    } elsif ($opt->compare && $done) {
        Analyze($analyze, $testDirName);
    } else {
        print "Running close annotation: $closeAnno\n";
        my $start = time;
        my @results = `$closeAnno`;
        my $duration = time - $start;
        print "$duration seconds to annotate.\n";
        open(my $oh, '>', "$testDirName/results.tbl") || die "Could not open results.tbl: $!";
        print $oh @results;
        close $oh; undef $oh;
        Analyze($analyze, $testDirName);
        print "Writing the completion flag.\n";
        open($oh, '>', "$testDirName/DONE") || die "Could not write completion file: $!";
        print $oh "$duration";
        close $oh; undef $oh;
    }
}

## Run the analysis for a single test.
sub Analyze {
    my ($analyze, $testDirName) = @_;
    print "Running the analysis:  $analyze\n";
    my @differences = `$analyze`;
    open(my $oh, '>', "$testDirName/compare.tbl") || die "Could not open compare.tbl: $!";
    print $oh @differences;
    close $oh; undef $oh;
}
