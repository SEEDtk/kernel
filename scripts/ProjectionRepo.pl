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
use Shrub;
use ScriptUtils;
use File::Copy::Recursive;
use Stats;
use BasicLocation;

=head1 Create Projection Repository

    ProjectionRepo.pl [ options ] outDir

This script creates the Projection Repository described in L<ProjectionChallenge>. The repository will be restricted to the current
set of well-behaved genomes in the Shrub. As currently constructed, it will only contain PEG features.

=head2 Parameters

The positional parameter is the name of the output directory. A C<GenomeData> and a C<SubsystemsData> directory will be created
underneath this named output directory.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item clear

If the output directory already exists, erase it before starting. The default is to fail if the output directory already exists.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir',
        Shrub::script_options(),
        ['clear|C', 'erase output directory before creating the repository']
        );
# Create the statistics object.
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Find the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (-d $outDir) {
    if ($opt->clear) {
        print "Erasing $outDir.\n";
        File::Copy::Recursive::pathempty($outDir);
    } else {
        die "Directory $outDir already exists.";
    }
} else {
    print "Creating $outDir.\n";
    create_dir($outDir);
}
# Locate the DNA repository.
my $dnaRepo = $shrub->DNArepo();
# Create the GenomeData directory.
print "Creating GenomeData.\n";
my $genomeRoot = "$outDir/GenomeData";
create_dir($genomeRoot);
# Get all the well-behaved genomes.
print "Loading genome list.\n";
my %genomes = map { $_->[0] => [$_->[1], $_->[2]] } $shrub->GetAll('Genome', 'Genome(well-behaved) = ?', [1], 'id name contig-file');
my $gCount = scalar keys %genomes;
print "$gCount genomes found.\n";
# Loop through the genomes.
my $gNum = 0;
for my $genome (sort keys %genomes) {
    ++$gNum;
    my ($gName, $contigFile) = @{$genomes{$genome}};
    print "Processing $genome $gName ($gNum of $gCount).\n";
    # Create the genome's directory.
    my $genomeDir = "$genomeRoot/$genome";
    create_dir($genomeDir);
    # Output the type file.
    create_list("$genomeDir/name", [$gName]);
    # Copy the contig file.
    File::Copy::Recursive::fcopy("$dnaRepo/$contigFile", "$genomeDir/contigs");
    $stats->Add(contigFilesCopied => 1);
    # Now we output the pegs. First, we open the files.
    open(my $ph, ">$genomeDir/peg-info") || die "Could not open peg-info file for $genome: $!";
    open(my $fh, ">$genomeDir/peg-trans") || die "Could not open peg-trans file for $genome: $!";
    print "Reading peg functions.\n";
    # Get all the pegs.
    my $q = $shrub->Get('Feature Feature2Function Function AND Feature Protein',
            'Feature(id) LIKE ? AND Feature2Function(security) = ?', ["fig|$genome.peg.%", 2],
            'Feature(id) Function(description) Protein(sequence)');
    print "Processing pegs.\n";
    while (my $pegData = $q->Fetch()) {
        my ($fid, $fun, $prot) = $pegData->Values([qw(Feature(id) Function(description) Protein(sequence))]);
        $stats->Add(pegsProcessed => 1);
        # Get the location of this feature.
        my @locs = $shrub->fid_locs($fid);
        my $locString = join(',', map { $_->String } @locs);
        # Output the two records for this peg.
        print $ph join("\t", $fid, $locString, $fun) . "\n";
        print $fh ">$fid\n$prot\n";
    }
    close $ph;
    close $fh;
    $stats->Add(genomesProcessed => 1);
}
# Now we need to process the subsystem.
print "Creating SubsystemsData.\n";
my $subsysRoot = "$outDir/SubsystemsData";
create_dir($subsysRoot);
# Get the subsystems and their roles.
print "Reading subsystem roles.\n";
my %subs;
my %subNames;
my $q = $shrub->Get('Subsystem Role', '', [], 'Subsystem(id) Subsystem(name) Role(description)');
while (my $subData = $q->Fetch()) {
    my ($sub, $subName, $role) = $subData->Values([qw(Subsystem(id) Subsystem(name) Role(description))]);
    push @{$subs{$sub}}, $role;
    $subNames{$sub} = $subName;
    $stats->Add(subRoleFound => 1);
}
my $subCount = scalar(keys %subs);
my $subNum = 0;
print "$subCount subsystems found.\n";
# Loop through the subsystems.
for my $sub (sort keys %subs) {
    ++$subNum;
    print "Processing $sub ($subNum of $subCount).\n";
    # This will track the genomes we've found.
    my %subGenomes;
    # Create the subsystem directory.
    my $subDir = "$subsysRoot/$subNames{$sub}";
    $subDir =~ tr/ /_/;
    create_dir($subDir);
    # Write the roles.
    create_list("$subDir/Roles", $subs{$sub});
    # Open the peg output file.
    open(my $oh, ">$subDir/PegsInSubsys") || die "Could not open peg output file for $sub: $!";
    # Find all the pegs in the subsystem.
    my @pegsInSubsys = $shrub->GetAll('Subsystem2Row Row2Cell Cell2Feature Feature2Function Function',
        'Subsystem2Row(from-link) = ? AND Feature2Function(security) = ?', [$sub, 2],
        'Feature2Function(from-link) Function(description)');
    # Loop through the pegs, writing the peg data.
    for my $pegData (@pegsInSubsys) {
        my ($peg, $function) = @$pegData;
        $stats->Add(subsysPegFound => 1);
        if ($peg =~ /^fig\|(\d+\.\d+)\.peg/ && $genomes{$1}) {
            my $genome = $1;
            $stats->Add(subsysPegKept => 1);
            $subGenomes{$genome}++;
            print $oh "$peg\t$function\n";
        }
    }
    my @genomesFound = sort keys %subGenomes;
    print scalar(@genomesFound) . " genomes use $sub.\n";
    create_list("$subDir/GenomesInSubsys", \@genomesFound);
}
print "All done.\n" . $stats->Show();

# Create a simple list file.
sub create_list {
    my ($fileName, $itemList) = @_;
    open(my $oh, ">$fileName") || die "Could not open $fileName: $!";
    $stats->Add(listFileCreated => 1);
    for my $item ($itemList) {
        print $oh "$item\n";
        $stats->Add(listItemOut => 1);
    }
    close $oh;
}

# Create a directory.
sub create_dir {
    my ($dirName) = @_;
    File::Copy::Recursive::pathmk($dirName) || die "Could not create directory $dirName: $!";
    $stats->Add(dirCreated => 1);
}
