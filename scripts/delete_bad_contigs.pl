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
use File::Copy::Recursive;
use RASTlib;
use GenomeTypeObject;

=head1 Create a New Genome Package By Deleting Bad Contigs From an Old One

    delete_bad_contigs.pl [ options ] packageDir package

This script operates on a genome package for which bad contigs have been identified using the L<find_bad_contigs.pl> script.
That script produces a tab-delimited file C<bad.contigs> in the package directory which contains bad contig IDs in the first
column and their length in the second column. The C<bin.fa> file will be read in and rewritten with those contigs removed.
This file will then be fed to RAST to create a GTO file. The two files will then be used to create a new genome package
and the ID of this package will be returned. The package ID will be written to the standard output; all other output will
be to STDERR.

This script will exit with a code of C<2> if the C<bad.contigs> file is not found.

=head2 Parameters

There are two positional parameters, the name and path of the genome packages directory, and the name of the individual
genome package to examine.

The command-line options are the following.

=over 4

=item user

User name for RAST access.

=item password

Password for RAST access.

=item sleep

Sleep interval in seconds while waiting for the job to complete. The default is C<60>.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir package', ScriptUtils::ih_options(),
        ["user|u=s", "user name for RAST access"],
        ["password|p=s", "password for RAST access"],
        ["sleep=i", "sleep interval for status polling", { default => 60 }],
        );
# Get the package directory and package.
my ($packageDir, $package) = @ARGV;
if (! $packageDir) {
    die "You must specified a package directory.";
} elsif (! -d $packageDir) {
    die "Invalid package directory $packageDir.";
} elsif (! $package) {
    die "You must specified a package.";
} elsif (! -d "$packageDir/$package") {
    die "Invalid package ID $package.";
} elsif (! -s "$packageDir/$package/bad.contigs") {
    print STDERR "find_bad_contigs was not run for $package.";
    exit(2);
}
my $pDir = "$packageDir/$package";
# Get a hash of the bad contigs.
open(my $ih, '<', "$pDir/bad.contigs") || die "Could not open bad.contigs: $!";
my %badContigs;
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^(\S+)\t/) {
        $badContigs{$1} = 1;
    }
}
close $ih; undef $ih;
# Figure out the names of the new files.
my $i = 1;
$i++ while (-f "$pDir/bin$i.fa");
my $newFa = "$pDir/bin$i.fa";
# Create the new FASTA file. We will also accumulate the contig triples in memory.
print STDERR "Spooling good contigs to $newFa.\n";
my @triples;
open($ih, '<', "$pDir/bin.fa") || die "Could not open bin.fa: $!";
open(my $oh, '>', $newFa) || die "Could not open $newFa: $!";
my ($keep, $contigID, $comment, @rows);
my ($contigs, $bases) = (0, 0);
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^>(\S+)(?:\s+(.+))?/) {
        # Here we are starting a new contig. Flush the old one and see if we want to keep the
        # new one.
        my ($newID, $newComment) = ($1, $2);
        if ($keep) {
            push @triples, [$contigID, $comment, join("", @rows)];
        }
        $contigID = $newID;
        $comment = $newComment || '';
        @rows = ();
        $keep = ($badContigs{$contigID} ? 0 : 1);
        if ($keep) {
            print $oh $line;
            $contigs++;
        }
    } elsif ($keep) {
        # Here we have a line of DNA data in a contig we're keeping.
        push @rows, $line;
        print $oh $line;
        chomp $line;
        $bases += length($line);
    }
}
if ($keep) {
    push @triples, [$contigID, $comment, join("", @rows)];
}
close $oh; undef $oh;
close $ih; undef $ih;
print STDERR "$contigs contigs found, $bases base pairs.\n";
# Get the taxonomy ID and name from the old GTO.
my ($taxonID, $name) = get_taxon($pDir);
$name = "$name (cleaned)";
# Create the GTO.
print STDERR "Submitting good contigs to RAST.\n";
my %options = ('sleep' => $opt->sleep, user => $opt->user, password => $opt->password);
my $gto = RASTlib::Annotate(\@triples, $taxonID, $name, %options);
# Output the GTO.
my $newID = $gto->{id};
print STDERR "Creating new GenomePackage $newID.\n";
my $newDir = "$packageDir/$newID";
if (-d $newDir) {
    die "Redundant genome ID returned $newID.";
} else {
    File::Copy::Recursive::pathmk($newDir) || die "Could not create $newDir: $!";
    File::Copy::Recursive::fmove($newFa, "$newDir/bin.fa") || die "Could not copy FASTA: $!";
    SeedUtils::write_encoded_object($gto, "$newDir/bin.gto");
    open($oh, '>', "$newDir/data.tbl") || die "Could not open new data.tbl: $!";
    print $oh "Genome Name\t$name\n";
    print $oh "Source Package\t$package\n";
    print $oh "Contigs\t$contigs\n";
    print $oh "Base pairs\t$bases\n";
    print "$newID\n";
}


sub get_taxon {
    my ($pDir) = @_;
    my $gto = GenomeTypeObject->create_from_file("$pDir/bin.gto");
    my $taxon = $gto->{ncbi_taxonomy_id};
    my $name = $gto->{scientific_name};
    return ($taxon, $name);
}