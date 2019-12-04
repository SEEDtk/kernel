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
working directory, then run L<Close_Anno.pl> and L<p3x-genome-function-file.pl>.

=head2 Parameters

The positional parameter is the name of the testing directory.

The command-line options are as follows:

=over 4

=item missing

Only run a test if it has not already been run.

=item compare

Only run the comparisons.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use GenomeTypeObject;
use IPC::Run3;
use File::Copy::Recursive;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('testDir',
        ['compare', 'rerun the comparisons'],
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
# Get the test genome files.
opendir(my $dh, $testDir) || die "Could not open $testDir: $!";
my @gFiles = sort grep { $_ =~ /^\d+\.\d+\.tbl$/ } readdir $dh;
print scalar(@gFiles) . " test files found.\n";
# Make the test directory current.
chdir "$testDir";
# This hash will map genome IDs to protein FASTA file names.
my %genomeFasta;
# This hash will map genome IDs to genetic codes.
my %genomeGC;
# Loop through the test genome files.
for my $gFile (@gFiles) {
    # Extract the genome ID.
    my ($genomeID) = ($gFile =~ /^(\d+\.\d+)/);
    print "Processing tests for $genomeID.\n";
    # Get the GTO for this genome and remove all the feature data.
    CreateSourceGto($genomeID);
    # Now read the neighbor file.
    my $nList = FindCloseGenomes($genomeID);
    my $nCount = scalar @$nList;
    if ($nCount < 2) {
        die "Neighbor file for $genomeID too small.";
    } else {
        print "$nCount neighbors found for $genomeID.\n";
        # Compute one-half the neighbor count.
        my $half = ($nCount + 1) >> 1;
        # Now we run the tests.
        RunTest($genomeID, $nList, 0, $half, "closest$half");
        RunTest($genomeID, $nList, $nCount - $half, $half, "furthest$half");
        RunTest($genomeID, $nList, 0, $nCount, "closest$nCount");
        if ($half > 2) {
            RunTest($genomeID, $nList, 0, 2, "closest2");
        }
        if ($half > 1) {
            RunTest($genomeID, $nList, $nCount - 1, 1, "furthest1");
        }
    }
}

sub CreateSourceGto {
    my ($genomeID) = @_;
    my $gtoName = "$genomeID.gto";
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

sub FindCloseGenomes {
    my ($genomeID) = @_;
    my @retVal;
    print "Reading neighbor file for $genomeID.\n";
    open(my $ih, '<', "$genomeID.tbl") || die "Could not open neighbor file for $genomeID: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /(\d+\.\d+)/) {
            my $neighborID = $1;
            if ($genomeFasta{$neighborID}) {
                print "$neighborID is a duplicate.\n";
            } else {
                print "Downloading neighbor $neighborID.\n";
                my $gto = $p3->gto_of($neighborID);
                if (! $gto) {
                    print "$neighborID not found-- skipping.\n";
                } else {
                    # Save the genetic code and create the FASTA.
                    $genomeGC{$neighborID} = $gto->{genetic_code} // 11;
                    print "Creating FASTA file for $neighborID.\n";
                    my $fastaFile = "$neighborID.fna";
                    $gto->write_protein_translations_to_file($fastaFile, 1);
                    $genomeFasta{$neighborID} = $fastaFile;
                    push @retVal, $neighborID;
                }
            }
        }
    }
    return \@retVal;
}

sub RunTest {
    my ($genomeID, $nList, $i0, $n, $testName) = @_;
    # Build the commands.
    my $testDirName = "$genomeID.$testName";
    my $closeAnno = "Close_Anno --log $testDirName/results.log --gto $testDirName/results.gto $genomeID.gto $testDirName/control.tbl $testDirName";
    my $analyze = "p3x-genome-function-file --log $testDirName/compare.log --input $testDirName/results.tbl $genomeID";
    # Check to see if we need to run this test.
    my $done = -f "$testDirName/DONE";
    if ($opt->missing && $done) {
        print "Skipping $testName for $genomeID:  already found.\n";
    } elsif ($opt->compare && $done) {
        Analyze($analyze, $testDirName);
    } else {
        # Create/clear the test directory.
        if (-d $testDirName) {
            print "Clearing $testDirName.\n";
            File::Copy::Recursive::pathempty($testDirName);
        } else {
            print "Creating $testDirName.\n";
            File::Copy::Recursive::pathmk($testDirName);
        }
        # Create the control file for this test.
        open(my $oh, '>', "$testDirName/control.tbl") || die "Could not create control file: $!";
        for (my $i = $i0; $i < $i0 + $n; $i++) {
            my $neighborID = $nList->[$i];
            print $oh "$genomeFasta{$neighborID}\t$genomeGC{$neighborID}\n";
        }
        close $oh; undef $oh;
        print "Running close annotation: $closeAnno\n";
        my $start = time;
        my @results = `$closeAnno`;
        print ((time - $start) . " seconds to annotate.\n");
        open($oh, '>', "$testDirName/results.tbl") || die "Could not open results.tbl: $!";
        print $oh @results;
        close $oh; undef $oh;
        Analyze($analyze, $testDirName);
        print "Writing the completion flag.\n";
        open($oh, '>', "$testDirName/DONE") || die "Could not write completion file: $!";
        print $oh "\n";
        close $oh; undef $oh;
    }
}

sub Analyze {
    my ($analyze, $testDirName) = @_;
    print "Running the analysis:  $analyze\n";
    my @differences = `$analyze`;
    open(my $oh, '>', "$testDirName/compare.tbl") || die "Could not open compare.tbl: $!";
    print $oh @differences;
    close $oh; undef $oh;
}
