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
use P3DataAPI;
use File::Copy::Recursive;


=head1 Create Genome Package from PATRIC Genome

    package_genome.pl [ options ] genomeID packageDir

This script will create a genome package from a PATRIC genome. Each genome package is a directory with the
same name as the genome ID and contains the following files

=over 4

=item bin.gto

A L<GenomeTypeObject> for the genome.

=item bin.fa

A FASTA file containing the genome.

=item data.tbl

A tab-delimited file of key/value pairs, with the following keys.

=over 8

=item Genome Name

Name of the genome

=item Source Database

PATRIC

=item Contigs

The number of contigs in the genome.

=item Base pairs

The number of DNA base pairs in the genome's contigs.

=back

=back

=head2 Parameters

The two positional parameters are the ID of the PATRIC genome and the full name of the directory to contain the output
genome packages. If the genome package already exists, it will not be overwritten.

The command-line options are the following.

=over 4

=item force

If specified, existing genome packages will be overwritten.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('genomeID packageDir',
        ['force', "overwrite existing packages"],
        );
# Get the directories.
my ($genomeID, $packageDir) = @ARGV;
if (! $genomeID) {
    die "No genome ID specified.";
} elsif (! ($genomeID =~ /\d+\.\d+/)) {
    die "Invalid genome ID $genomeID.";
}
if (! $packageDir) {
    die "No package directory specified.";
} elsif (! -d $packageDir) {
    die "Invalid package directory $packageDir";
}
# Connect to the database. Note that we must do an environment hack.
$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME} = 0;
my $d = P3DataAPI->new();
# Get the GTO.
my $gto = $d->gto_of($genomeID);
if (! $gto) {
    die "Genome $genomeID not found.";
} else {
    # Get the name of the genome.
    my $name = $gto->{scientific_name};
    # Compute the output directory name.
    my $genomeDir = "$packageDir/$genomeID";
    if (-d $genomeDir && ! $opt->force) {
        print "Package already exists for $genomeID.\n";
    } else {
        # Here the bin is new or we are forcing.
        if (! -d $genomeDir) {
            print "Creating $genomeDir.\n";
            mkdir $genomeDir;
        } else {
            print "Replacing $genomeDir.\n";
            File::Copy::Recursive::pathempty($genomeDir);
        }
        # Loop through the contigs, creating the FASTA file and counting contigs and DNA.
        open(my $fh, '>', "$genomeDir/bin.fa") || die "Could not open FASTA output: $!";
        my $contigs = $gto->{contigs};
        my $contigCount = scalar @$contigs;
        my $dnaLen = 0;
        for my $contig (@$contigs) {
            my $dna = $contig->{dna};
            $dnaLen += length $dna;
            print $fh ">$contig->{id}\n$dna\n";
        }
        close $fh;
        # Output the GTO file.
        $gto->destroy_to_file("$genomeDir/bin.gto");
        # We will compute the data table values in here.
        my %data;
        $data{'Genome Name'} = $name;
        $data{'Source Database'} = 'PATRIC';
        $data{'Contigs'} = $contigCount;
        $data{'Base pairs'} = $dnaLen;
        # Find the closest reference genome.
        # Write the data file.
        open(my $oh, '>', "$genomeDir/data.tbl") || die "Could not open data table file: $!";
        for my $key (sort keys %data) {
            print $oh "$key\t$data{$key}\n";
        }
        print "All done.\n";
    }
}