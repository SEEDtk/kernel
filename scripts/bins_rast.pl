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
use Bin;
use Bin::Analyze;
use Loader;
use Shrub;

=head1 Submit Bins Through RAST

    bins_rast.pl [ options ] workDir

This script submits a set of bins through RAST to create proposed genomes from them. The C<bins.json> file is
read to determine the bin contents, and then the C<sample.fasta> file is processed to create the FASTA files.
The C<ref.genomes.scores.tbl> is used to compute the names of the genomes (if present).

=head2 Parameters

The single positional parameter is the working directory used to access the bin data.

The command-line options are the following.

=over 4

=item user

User name for RAST access. If omitted, the default is taken from the RASTUSER environment variable.

=item password

Password for RAST access. If omitted, the default is taken from the RASTPASS environment variable.

=item sleep

Sleep interval in seconds while waiting for the job to complete. The default is C<60>.

=item partial

Only call RAST for bins for which a GTO does not yet exist. Use this to resume after a failure.

=back

=head2 Output

A JSON-format L<GenomeTypeObject> will be produced for each bin, with the name C<bin>I<X>C<.gto>, with
I<X> being the bin number.

A new copy of the bin file with universal role information embedded will be output to C<bins.rast.json>.

The analysis report will be updated in C<bins.report.txt> with the universal role data.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir', Shrub::script_options(),
        ["user|u=s", "user name for RAST access"],
        ["password|p=s", "password for RAST access"],
        ["sleep=i", "sleep interval for status polling", { default => 60 }],
        ['partial', 'only process new bins'],
        );
# Verify the work directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No work directory specified.";
} elsif (! -d $workDir) {
    die "Invalid work directory $workDir specified.";
} else {
    my $binJsonFile = "$workDir/bins.json";
    my $contigFastaFile = "$workDir/sample.fasta";
    if (! -s $binJsonFile || ! -s $contigFastaFile) {
        die "$workDir does not appear to contain completed bins.";
    }
    # Connect to the database.
    my $shrub = Shrub->new_for_script($opt);
    # Create a hash of the universal roles.
    my $uniRoleH = $shrub->GetUniRoles();
    my $totUnis = scalar keys %$uniRoleH;
    # Create an analyzer object.
    my $analyzer = Bin::Analyze->new($shrub, totUnis => $totUnis, minUnis => (0.8 * $totUnis));
    # Get the loader object.
    my $loader = Loader->new();
    my $stats = $loader->stats;
    # Read in the bins.
    print "Reading bins from $binJsonFile.\n";
    my $binList = Bin::ReadBins($binJsonFile);
    # Create the RAST option hash.
    my %rastOpts = (user => $opt->user, password => $opt->password, 'sleep' => $opt->sleep);
    # Loop through the bins, processing them one at a time. For each bin, we read the whole sample
    # file to get the contigs. Only one bin's worth of contigs is kept in memory at a time, at the cost of
    # a lot of extra reading.
    my $binNum = 0;
    for my $bin (@$binList) {
        my $taxonID = $bin->taxonID;
        my $name = $bin->name;
        $stats->Add(bins => 1);
        $binNum++;
        print "Processing bin $binNum - $taxonID: $name.\n";
        my $binFile = "$workDir/bin$binNum.gto";
        my $gto;
        if ($opt->partial && -s $binFile) {
            print "Reading bin GTO from $binFile.\n";
            $stats->Add(binsSkipped => 1);
            $gto = SeedUtils::read_encoded_object($binFile);
        } else {
            my %contigs = map { $_ => 1 } $bin->contigs;
            # Now we read the sample file and keep the contig triples.
            my @triples;
            my $ih = $loader->OpenFasta(sampleContigs => $contigFastaFile);
            open(my $oh, ">$workDir/bin$binNum.fa") || die "Could not open FASTA output file for bin $binNum.";
            my $triple = $loader->GetLine(sampleContigs => $ih);
            while (defined $triple) {
                my $contigID = $triple->[0];
                if ($contigs{$contigID}) {
                    $stats->Add(contigsKept => 1);
                    push @triples, $triple;
                    print $oh ">$triple->[0] $triple->[1]\n$triple->[2]\n";
                }
                $triple = $loader->GetLine(sampleContigs => $ih);
            }
            my $contigCount = scalar @triples;
            print "Submitting $binNum to RAST: $contigCount contigs.\n";
            $gto = RASTlib::Annotate(\@triples, $taxonID, $name, %rastOpts);
            print "Spooling genome to $workDir.\n";
            SeedUtils::write_encoded_object($gto, "$workDir/bin$binNum.gto");
        }
        print "Searching for universal proteins.\n";
        # Clear the bin's current universal protein list.
        $bin->replace_prots();
        # Search the genome for universal roles.
        my $flist = $gto->{features};
        for my $feature (@$flist) {
            my $fun = $feature->{function};
            if ($fun) {
                my $funID = $shrub->desc_to_function($fun);
                if ($funID && $uniRoleH->{$funID}) {
                    $bin->incr_prot($funID, 1);
                }
            }
        }
        my $protH = $bin->uniProts;
        my $protCount = scalar keys %$protH;
        print "$protCount universal proteins found.\n";
    }
    # Output the new bins.
    my $binOutputFile = "$workDir/bins.rast.json";
    print "Spooling bins to $binOutputFile.\n";
    my @sorted = sort { Bin::cmp($a, $b) } @$binList;
    open(my $oh, ">$binOutputFile") || die "Could not open bins.rast.json output file: $!";
    for my $bin (@sorted) {
        $bin->Write($oh);
    }
    close $oh;
    # Read the reference genome file. We need this for the report.
    my $refScoreFile = "$workDir/ref.genomes.scores.tbl";
    if (-s $refScoreFile) {
        my %genomes;
        print "Reading reference genomes from $refScoreFile.\n";
        my $ih = $loader->OpenFile(refGenomes => $refScoreFile);
        while (my $refFields = $loader->GetLine(refGenomes => $ih)) {
            my ($contig, $genome, $score, $name) = @$refFields;
            $genomes{$genome} = $name;
        }
        $analyzer->SetGenomes(\%genomes);
    }
    # Write the report.
    open(my $rh, ">$workDir/bins.report.txt") || die "Could not open report file: $!";
    $analyzer->BinReport($rh, $uniRoleH, \@sorted);
    close $rh;
}
print "All done.\n";
