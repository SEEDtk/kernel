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
use RASTlib;
use TaxCheck;
use File::Copy::Recursive;
use P3DataAPI;
use Contigs;
use SeedUtils;
use GEO;
use EvalCon;

=head1 RAST Scramble for METABAT Test

    rast_scramble.pl [ options ] outDir

This script runs a bunch of RAST jobs in parallel.  It takes as input a 2-column file.  The first column is a FASTA file and the second column is a
taxonomy ID.  If the second column is "66666", the taxonomy ID must be estimated.

The last part of the path in the FASTA file name is used to compute a subdirectory in the output directory.  The true file name is used to compute a
file name in that subdirectory.  A RAST job is then run to annotate the FASTA file.  We will start a fixed number of jobs, and then retrieve them all,
then repeat the process.

=head2 Parameters

The positional parameter is the name of the output directory.  The command-line options are those L<ScriptUtils/ih_options> plus the following:

=over 4

=item clear

Erase the output directory before processing.

=item processes

Number of processes to run in parallel.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir', ScriptUtils::ih_options(),
        ['clear', 'erase output before starting'],
        ['processes=i', 'number of processes to run in each batch']
        );
# [domain, gc] information will be cached in here.
my %taxons;
# Create a statistics object.
my $stats = Stats->new();
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Get the taxonomy checker.
my $p3 = P3DataAPI->new();
print STDERR "Connecting to taxonomy checker.\n";
my $taxChecker = TaxCheck->new($p3, debug => 1);
# Load the role hashes.
print STDERR "Loading role hashes.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
# Insure we have an output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print STDERR "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory: $!";
} elsif ($opt->clear) {
    print STDERR "Erasing output directory $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase output directory: $!";
} else {
    print STDERR "Writing output to $outDir.\n"
}
# Get the RAST authorization header.
my $auth_header = RASTlib::auth_header();
# This hash will count the good and bad genomes in each subdirectory as name => [good, bad].
my %counts;
# Loop through the input.
while (! eof $ih) {
    # Get the next batch of FASTA files.
    my @couplets = ScriptUtils::get_couplets($ih, 1, $opt->processes);
    $stats->Add(batchIn => 1);
    # This list will contain the job ID, checksum, taxon ID, subDir, and output file name for each genome.  If the
    # output file already exists, only the output file name will be defined.
    my @batch;
    for my $couplet (@couplets) {
        my ($file, $data) = @$couplet;
        my $taxID = $data->[1];
        $stats->Add(fileIn => 1);
        # Compute the name and subdirectory.
        my @parts = split '/', $file;
        my $name = pop @parts;
        my $subDir = pop @parts;
        my $outSubDir = "$outDir/$subDir";
        if (! -d $outSubDir) {
            File::Copy::Recursive::pathmk($outSubDir) || die "Could not create $outSubDir: $!";
        }
        my $outFile = "$outSubDir/$name.gto";
        if (-s $outFile) {
            print STDERR "Output GTO $outFile already exists.\n";
            push @batch, [undef, undef, undef, $subDir, $outFile];
            $stats->Add(outputAlreadyPresent => 1);
        } else {
            print STDERR "Annotating $file into $outFile.\n";
            # Insure we have a good taxonomy ID.
            if ($taxID eq '66666') {
                ($taxID) = $taxChecker->Compute($file, 'species');
                $stats->Add(taxIdComputed => 1);
                $taxID //= '66666';
            }
            my ($domain, $gc);
            if ($taxons{$taxID}) {
                ($domain, $gc) = @{$taxons{$taxID}};
            } else {
                $stats->Add(taxDataComputed => 1);
                my ($result) = $p3->query('taxonomy', ['select' => 'division', 'genetic_code'],
                        ['eq' => 'taxon_id', $taxID]);
                if (! $result) {
                    print STDERR "WARNING: taxonomic ID $taxID not found in PATRIC.\n";
                    ($domain, $gc) = ('Bacteria', 11);
                } else {
                    ($domain, $gc) = ($result->{division}, $result->{genetic_code});
                    if ($domain ne 'Archaea' && $domain ne 'Bacteria') {
                        print STDERR "WARNING: taxonomic ID $taxID has nonstandard domain $domain.\n";
                        $domain = 'Bacteria';
                        $stats->Add(badDomain => 1);
                    }
                }
                $taxons{$taxID} = [$domain, $gc];
            }
            # Read in the FASTA file.
            my $contigs = Contigs->new($file);
            # Submit the annotation.
            print STDERR "Submitting to RAST with taxonomic ID $taxID in domain $domain with gc = $gc.\n";
            my ($jobID, $checksum) = RASTlib::Submit([$contigs->tuples], $taxID, $name, domain => $domain,
                    geneticCode => $gc, noIndex => 1, robust => 1);
            if ($jobID) {
                push @batch, [$jobID, $checksum, $taxID, $subDir, $outFile];
                $stats->Add(jobSubmitted => 1);
            } else {
                $stats->Add(submitFailed => 1);
            }
        }
    }
    # Now all the jobs are submitted.  Get them back.
    while (scalar @batch) {
        # We pass through the whole batch, returning the incomplete jobs to this list for the
        # next pass.
        my @remaining;
        for my $job (@batch) {
            my ($jobID, $checksum, $taxID, $subDir, $outFile) = @$job;
            if (! $jobID) {
                # Here we have to read the GTO from disk.
                my $gto = GenomeTypeObject->create_from_file($outFile);
                qual_check($subDir, $gto);
                $stats->Add(gtoRead => 1);
            } else {
                # Here the GTO has to be retrieved from RAST.
                my $completed = RASTlib::check($jobID, header => $auth_header, robust => 1);
                if ($completed < 0) {
                    print STDERR "WARNING: job $jobID for $outFile failed.\n";
                    $stats->Add(rastError => 1);
                } elsif ($completed == 0) {
                    push @remaining, $job;
                } else {
                    # Here the job is done.
                    my $gto = RASTlib::retrieve($jobID, $checksum, $auth_header, $taxID);
                    print STDERR "Writing output for $jobID to $outFile.\n";
                    qual_check($subDir, $gto);
                    $gto->destroy_to_file($outFile);
                    $stats->Add(gtoWritten => 1);
                }
            }
        }
        # Prepare for the next run through.
        @batch = @remaining;
        # Only check every 10 seconds.
        if (scalar @batch) {
            sleep 10;
        }
    }
    # End of batch, process the next one.
}
# All done. Write the quality report.
print STDERR "Writing counts to $outDir/counts.tbl\n";
open(my $oh, '>', "$outDir/counts.tbl") || die "Could not open counts.tbl: $!";
print $oh "Sample\tGood\tBad\n";
for my $subDir (sort keys %counts) {
    print $oh "$subDir\t$counts{$subDir}[0]\t$counts{$subDir}[1]\n";
}
close $oh;
print "All done.\n" . $stats->Show();

sub qual_check {
    my ($subDir, $gto) = @_;
    my $geo = GEO->CreateFromGto($gto, roleHashes => [$nMap, $cMap], detail => 0, logH => \*STDERR, p3 => $p3);
    my $counters = $counts{$subDir};
    if (! $counters) {
        $counters = [0, 0];
        $counts{$subDir} = $counters;
    }
    if ($geo->is_good) {
        $stats->Add(goodGenome => 1);
        $counters->[0]++;
    } else {
        $stats->Add(badGenome => 1);
        $counters->[1]++;
    }
}
