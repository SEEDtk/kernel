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

=head1 Generate Role-Triple Files

    p3-role-file-set.pl [options] outDir role1 role2 ... roleN

This script handles the incredibly complex tax of generating role-triple files.  These files are generated for good
genomes in BV-BRC for one or more roles.  The process works in two passes for each role.  During the first pass, we
pull in the protein strings and remove bad genomes or genomes with more than one copy of the protein.  In the second
pass we remove proteins whose length differs too much from the mean.

Finally, we read the files back in and eliminate genomes that are not in all the files.

The initial output files will have the name I<roleName>C<.raw.tbl>.  The final output files will have the name
I<roleName>C<.rep.tbl>.  Each output file will have as columns the genome ID, the quality score, and the protein
sequence.

=head2 Parameters

The positional parameters are the name of the output directory and then the IDs of the roles.  The roles must be
present in the standard C<roles.in.subsystems> file or in an equivalent file specified as input.

The command-line options are as follows.

=over 4

=item roleFile

The C<roles.in.subsystems> file containing the role definitions.  The default is taken from the SEEDtk global data
directory.

=item goodFile

The C<patric.good.tbl> file containing the good genomes and taxonomic lineages.  The default is taken from the SEEDtk
data directory.

=item preFiltered

The ID of a role that has already been length-filtered in the good genome list.  The default is C<PhenTrnaSyntAlph>.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use SeedUtils;
use Stats;
use File::Copy::Recursive;
use Statistics::Descriptive;
use RoleParse;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir role1 role2 ... roleN',
        ['roleFile|roles|r=s', 'file containing role definitions', { default => "$FIG_Config::p3data/roles.in.subsystems" }],
        ['goodFile|good|g=s', 'file containing good genome definitions', { default => "$FIG_Config::data/patric.good.tbl"}],
        ['preFiltered|pre|p=s', 'prefiltered role', { default => 'PhenTrnaSyntAlph' }],
    );
my $stats = Stats->new();
# Connect to BV-BRC.
my $p3 = P3DataAPI->new();
# Get the options.
my $roleFile = $opt->rolefile;
my $goodFile = $opt->goodfile;
my $preFRole = $opt->prefiltered;
# Get the good genomes.  For each one, we need to know its domain (A or B), its score, and the number of roles for
# which it was written.
my (%gDomain, %gScore, %gRoles);
print "Reading $goodFile for domains and scores.\n";
open(my $ih, '<', $goodFile) || die "Could not open good-genome file: $!";
my ($headers, $cols) = P3Utils::find_headers($ih, goodGenomes => 'genome_id', 'taxon_lineage_ids', 'Score');
while (! eof $ih) {
    my ($genome, $lineage, $score) = P3Utils::get_cols($ih, $cols);
    $stats->Add(goodGenome => 1);
    if ($lineage =~ /::2157::/) {
        $stats->Add(goodArchaea => 1);
        $gDomain{$genome} = 'A';
    } elsif ($lineage =~ /::2::/) {
        $stats->Add(goodBacteria => 1);
        $gDomain{$genome} = 'B';
    }
    $gScore{$genome} = $score;
    $gRoles{$genome} = 0;
}
close $ih; undef $ih;
# Get the positional parameters.
my ($outDir, @repRoles) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive($outDir) || die "Could not create directory: $!";
}
if (! @repRoles) {
    die "No role IDs specified.";
}
# Now create the role hashes.  For each role we need its name and checksum.
print "Loading role definitions from $roleFile.\n";
my %rName = map { $_ => '' } @repRoles;
my %rCheck;
open($ih, '<', $roleFile) || die "Could not open roleFile: $!";
while (! eof $ih) {
    my ($roleID, $checksum, $name) = P3Utils::get_fields($ih);
    if (exists $rName{$roleID}) {
        $rName{$roleID} = $name;
        $rCheck{$roleID} = $checksum;
    }
}
close $ih; undef $ih;
my @notFound = grep { ! $rName{$_} } @repRoles;
if (@notFound) {
    die "Missing roles in definition file: " . join(", ", @notFound) . ".\n";
}
# We have all our pieces in place.  The next step is to loop through the creation of the raw triple files.
for my $repRole (@repRoles) {
    my $rawFile = "$outDir/$repRole.raw.tbl";
    print "Reading proteins for $repRole.\n";
    # Initialize the stats computers.
    my %computer = map { $_ => Statistics::Descriptive::Sparse->new() } qw(A B);
    # This hash will hold the protein strings.
    my %gProt;
    # Get and clean the role name.
    my $protName = $rName{$repRole};
    die "$repRole is an invalid role ID." if ! $protName;
    $protName =~ s/\([^)]+\)//g;
    $protName =~ s/\s+$//;
    $protName =~ s/\s+/ /g;
    # Get the role checksum.
    my $roleCheck = $rCheck{$repRole};
    # Get all the protein instances.
    my $proteins = P3Utils::get_data($p3, feature => [['eq', 'product', $protName]], ['genome_id', 'product', 'aa_sequence']);
    my $totFound = scalar @$proteins;
    print "$totFound proteins found.\n";
    for my $protein (@$proteins) {
        my ($genome, $function, $seq) = @$protein;
        if (! $genome) {
            $stats->Add(missingGenome => 1);
        } elsif (! $seq) {
            $stats->Add(missingSequence => 1);
        } elsif (! $gDomain{$genome}) {
            # Not one of the genomes in the good-genome list.
            $stats->Add(badGenome => 1);
        } elsif (RoleParse::Checksum($function) ne $roleCheck) {
            # Not our role, but has enough words in common that it is a SOLR match.
            $stats->Add(funnyRole => 1);
        } elsif (exists $gProt{$genome}) {
            $stats->Add(dupRole => 1);
            # Insure we don't use this sequence for anything.  The protein is most likely broken and would throw off
            # the mean and standard deviation.
            $gProt{$genome} = '';
        } else {
            $gProt{$genome} = $seq;
            $stats->Add(seqStored => 1);
        }
    }
    # Release the memory for all the proteins.
    undef $proteins;
    # We have all the protein sequences.  The next step is to compute the length limits.
    my %limits;
    if ($repRole eq $preFRole) {
        # This is a pre-filtered role, so we just put in ridiculous limits.
        %limits = (A => [0, 1000000], B => [0, 1000000]);
    } else {
        # Here we have to scan to compute the mean and standard deviation.
        print "Computing length limits for $repRole.\n";
        # Allocate the statistical objects.
        %limits = map { $_ => Statistics::Descriptive::Sparse->new() } qw(A B);
        # Loop through the sequences.
        for my $genome (keys %gProt) {
            if ($gProt{$genome}) {
                $limits{$gDomain{$genome}}->add_data(length $gProt{$genome});
            }
        }
        # Compute the min and max using the mean and standard deviation.
        for my $d (keys %limits) {
            my $m = $limits{$d}->mean;
            my $s = $limits{$d}->standard_deviation * 3;
            print "Limit for [$d] $repRole length is $m +- $s.\n";
            $limits{$d} = [$m - $s, $m + $s];
        }
    }
    # Now we make a third pass to output the raw proteins filtering on length.
    # Here we also count the genome output in %gRoles;
    print "Writing proteins to $rawFile.\n";
    open(my $oh, '>', $rawFile) || die "Could not open raw output for $repRole: $!";
    P3Utils::print_cols(['genome_id', 'score', 'aa_sequence'], oh => $oh);
    for my $genome (sort keys %gProt) {
        my $seq = $gProt{$genome};
        my $len = length $seq;
        my $limits = $limits{$gDomain{$genome}};
        if (! $seq) {
            $stats->Add("$repRole-tooMany" => 1);
        } elsif ($len < $limits->[0]) {
            $stats->Add("$repRole-tooShort" => 1);
        } elsif ($len > $limits->[1]) {
            $stats->Add("$repRole-tooLong" => 1);
        } else {
            $stats->Add("$repRole-rawOut" => 1);
            P3Utils::print_cols([$genome, $gScore{$genome}, $seq], oh => $oh);
            $gRoles{$genome}++;
        }
    }
    close $oh; undef $oh;
}
# Now we go through the raw files to create the final files.  Only genomes whose counts equal the number of roles
# will be output.
my $target = scalar @repRoles;
print "Number of roles is $target.\n";
for my $repRole (@repRoles) {
    print "Finalizing $repRole.\n";
    open(my $ih, '<', "$outDir/$repRole.raw.tbl") || die "Could not open raw input for $repRole: $!";
    open(my $oh, '>', "$outDir/$repRole.rep.tbl") || die "Could not open final output for $repRole: $!";
    # Echo the header;
    my $line = <$ih>;
    print $oh $line;
    while (! eof $ih) {
        # Get the next line in the raw file and pull off the genome ID.
        $line = <$ih>;
        my ($genome) = split /\t/, $line;
        if ($gRoles{$genome} < $target) {
            $stats->Add("$repRole-rejected" => 1);
        } else {
            print $oh $line;
            $stats->Add("$repRole-final" => 1);
        }
    }
}
print "All done.\n" . $stats->Show();