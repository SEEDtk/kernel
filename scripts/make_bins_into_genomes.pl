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
use RASTlib;
use File::Copy::Recursive;
use Bin::Blast;
use Loader;
use SeedUtils;

=head1 Convert Bins to Genome Packages

    make_bins_into_genomes.pl [ options ] inDir packageDir sampleName

This script will convert the bins produced by L<assemble_bins.pl> into genome packages in the specified
target directory.

Each bin is represented by a file in the input directory with the name I<XXXXXX.X>C<.fasta>, where I<XXXXXX.X>
is the ID of the nearest representative genome. We will blast that genome's PhenTrnaSyntAlph protein against
the FASTA file to locate the gene, and then blast the gene against the PhenTrnaSyntAlph BLAST database in order
to find the closest genome. That genome's taxonomic ID will be used to name the genome and then invoke RAST to
annotate it. The resulting GTO will be used to create a genome package in the output package directory.

=head2 Parameters

The positional parameters are the name of the input directory containing the bins, the name of the output directory
for the genome packages, and the name of the sample from which the bins were generated. This last is used in the
metadata for the genome package.

The command-line options are those in L<Shrub::script_options> plus the following.

=over 4

=item user

The user name for logging into the PATRIC RAST. The default is taken from the C<RASTUSER> environment variable.

=item password

The password for logging into the PATRIC RAST. The default is taken from the C<RASTPASS> environment variable.

=item force

If specified, all bins will be processed. Otherwise, only bins that have not had packages created in the output
directory will be processed.

=item seedfasta

The name of the BLAST database for the seed protein in the various PATRIC genomes. The default is
C<PhenTrnaSyntAlph.fa> in the global data directory.

=item seedprot

The function ID of the seed protein. The default is C<PhenTrnaSyntAlph>. This must match the BLAST database
in the B<seedfasta> parameter.

=item maxE

The maximum acceptable E-value for BLAST results. The default is C<1e-20>.

=item refMaxE

The maximum acceptable E-value for reference genome BLAST results. The default is C<1e-10>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=back

=head2 Output

For each bin I<XXXXXX.X>C<.fasta>, a L<GenomeTypeObject> file will be created named I<XXXXXX.X><.gto>. Unlike
most GTO files, the file name is not taken from the described genome's ID, but from the reference genome ID.
A temporary file named C<prot.fa> will be created to hold the protein sequence for BLASTing.

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir packageDir sampleName',
        Shrub::script_options(),
        ["user|u=s", "user name for RAST access", { default => $ENV{RASTUSER} }],
        ["password|p=s", "password for RAST access", { default => $ENV{RASTPASS} }],
        ["force", "re-process previously processed bins"],
        ['seedfasta|F=s',  'BLAST database (or FASTA file) of seed protein in all genomes',
                           { default => "$FIG_Config::global/PhenTrnaSyntAlph.fa"}],
        ['maxE|e=f',       'maximum acceptable e-value for blast hits', { default => 1e-20 }],
        ['refMaxE=f',      'maximum acceptable e-value for reference genome blast hits', { default => 1e-10 }],
        ['gap|g=i',        'maximum permissible gap between blast hits for merging', { default => 600 }],
        ['minlen|l=f',     'minimum fraction of the protein that must match in a blast hit', { default => 0.5 }],
        ['seedprot=s',     'function ID of seed protein', { default => 'PhenTrnaSyntAlph'}],
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get a loader utility object.
my $loader = Loader->new();
# Extract the statistics object.
my $stats = $loader->stats;
# Get the BLAST options.
my $maxE = $opt->maxe;
my $gap = $opt->gap;
my $minlen = $opt->minlen;
my $refMaxE = $opt->refmaxe;
# Get the positional parameters.
my ($inDir, $packageDir, $sampleName) = @ARGV;
# Verify the directories.
if (! $inDir) {
    die "No input directory specified."
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
}
if (! $packageDir) {
    die "No output directory specified.";
} elsif (-f $packageDir) {
    die "Invalid output directory $packageDir specified.";
} elsif (! -d $packageDir) {
    File::Copy::Recursive::pathmk($packageDir) || die "Could not create output directory $packageDir: $!";
}
if (! $sampleName) {
    die "No sample name specified.";
}
# Compute the name of the protein FASTA file.
my $protFile = "$inDir/prot.fasta";
# Get the seed protein's functional ID and its BLAST database.
my $protBlastDb = $opt->seedfasta;
my $protID = $opt->seedprot;
# This counter will be used to assign a bin number.
my $binNum = 1;
# Get the input bins.
opendir(my $dh, $inDir) || die "Could not open input directory $inDir: $!";
my @inFiles = grep { $_ =~ /^\d+\.\d+\.fasta$/ } readdir $dh;
closedir $dh;
print scalar(@inFiles) . " FASTA files found in $inDir.\n";
# Loop through the files.
for my $inFile (sort @inFiles) {
    # Initialize the forcing flag.
    my $force = $opt->force;
    # Extract the reference genome ID. Because of the grep filter above we know this will work.
    my ($refGenomeID) = ($inFile =~ /^(\d+\.\d+)/);
    # Get its name from the database.
    my ($refName) = $shrub->GetFlat('Genome', 'Genome(id) = ?', [$refGenomeID], 'name');
    $refName //= 'Unknown genome $refGenomeID';
    print "Processing bin for $refGenomeID: $refName.\n";
    # Compute the file names.
    my $contigFile = "$inDir/$inFile";
    my $gtoFile = "$inDir/$refGenomeID.gto";
    # The actual GenomeTypeObject will go in here.
    my $gto;
    # Do we need to do this genome?
    if (-s $gtoFile && ! $force)  {
        $stats->Add(gtoAlreadyFound => 1);
        print "Genome $gtoFile already created for $refGenomeID bin.\n";
        # Read the GTO into memory.
        $gto = GenomeTypeObject->create_from_file($gtoFile);
    } else {
        # Get the seed protein in the reference genome.
        print "Retrieving seed protein from $refGenomeID.\n";
        my ($protInfo) = $shrub->GetAll('Function2Feature Feature2Protein Protein',
                'Function2Feature(from-link) = ? AND Function2Feature(to-link) LIKE ? AND Function2Feature(security) = ?',
                [$protID, "fig|$refGenomeID.peg.%", 1], 'Feature2Protein(from-link) Protein(sequence)');
        die "Seed protein not found in $refGenomeID." if ! $protInfo;
        # Create a FASTA file of the SEED protein.
        open(my $oh, ">$protFile") || die "Could not open protein output file: $!";
        print $oh ">$protInfo->[0]\n$protInfo->[1]\n";
        close $oh;
        # Load the bin's contigs into a BLAST object.
        my $blaster = Bin::Blast->new($shrub, $inDir, "$inDir/$inFile", maxE => $maxE, gap => $gap, minlen => $minlen);
        # Blast to find the seed protein.
        my $hitHash = $blaster->FindProtein($protFile);
        # Locate the longest hit.
        my ($hitCount, $bestLen, $bestHit) = (0, 0);
        for my $contig (keys %$hitHash) {
            my $hit = $hitHash->{$contig};
            if ($hit->Length >= $bestLen) {
                $bestHit = $hit;
                $bestLen = $hit->Length;
            }
            $hitCount++;
            $stats->Add(seedProteinHit => 1);
        }
        if (! $bestHit) {
            print "No seed protein found. Bin skipped.\n";
            $stats->Add(noSeedProtein => 1);
        } else {
            my $contig = $bestHit->Contig;
            print "Best of $hitCount hits is in contig $contig.\n";
            # Get the protein's DNA.
            my $seqHash = $loader->GetDNA({ $contig => $bestHit }, $contigFile);
            # Blast it to compute the closest genome.
            my $contigHash = $blaster->MatchProteins($seqHash, $protID, 1, $refMaxE, db => $protBlastDb, type => 'dna');
            my $closeGenomeList = $contigHash->{$contig};
            if (! $closeGenomeList || ! @$closeGenomeList) {
                print "No close genome found. Cannot classify bin $refGenomeID.\n";
                $stats->Add(noCloseGenome => 1);
            } else {
                # Pull out the genome ID of the best genome.
                my $closeGenome = $closeGenomeList->[0][0];
                # Get the taxonomic grouping.
                my ($taxID) = ($closeGenome =~ /^(\d+)/);
                print "Closest genome is $closeGenome.\n";
                my ($taxName) = $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(id) = ?', [$taxID], 'scientific-name');
                if (! $taxName) {
                    # The database does not have this taxonomic ID, so use the name of the closest genome.
                    $stats->Add(unknownTaxID => 1);
                    $taxName = $closeGenomeList->[0][2];
                }
                # Compute the genome name.
                my $genomeName = "$taxName from biosample $sampleName";
                # Load the contigs into memory.
                print "Loading contigs.\n";
                my @contigs;
                my $gh = $loader->OpenFasta(binContig => $contigFile);
                while (my $tuple = $loader->GetLine(binContig => $gh)) {
                    push @contigs, $tuple;
                }
                # Calculate the GTO.
                print "Invoking RAST service.\n";
                $gto = RASTlib::Annotate(\@contigs, $closeGenome, $genomeName, user => $opt->user, password => $opt->password);
                # Write the GTO to disk.
                SeedUtils::write_encoded_object($gto, $gtoFile);
                # Denote we are forcing creation of a package.
                $force = 1;
            }
        }
    }
    # Now we check to make sure we have a GTO from which we can create a package. Note that it will
    # be unblessed at this point, so we can use members but not methods.
    if ($gto) {
        # We have a GTO. Check to see if we need to create a package.
        my $genomeID = $gto->{id};
        print "Bin for $refGenomeID produced $genomeID.\n";
            my $genomePackage = "$packageDir/$genomeID";
        if (! $force && -s "$genomePackage/bin.gto") {
            print "Package for $genomeID already exists.\n";
            $stats->Add(packageSkipped => 1);
        } else {
            # Here we want to create a package. Get the key metadata from the GTO.
            my %data;
            $data{'Genome Name'} = $gto->{scientific_name};
            $data{'Sample Name'} = $sampleName;
            $data{'Bin Number'} = $binNum;
            $data{'Ref Genome'} = $refGenomeID;
            $data{'Ref Name'} = $refName;
            # Count the contigs and the DNA.
            my $contigs = $gto->{contigs};
            my $contigCount = scalar @$contigs;
            my $dnaLen = 0;
            for my $contig (@$contigs) {
                my $dna = $contig->{dna};
                $dnaLen += length $dna;
            }
            $data{'Contigs'} = $contigCount;
            $data{'Base pairs'} = $dnaLen;
            # Insure we have the output directory.
            if (! -d $genomePackage) {
                $stats->Add(packageAdded => 1);
                print "Creating $genomePackage directory.\n";
                File::Copy::Recursive::pathmk($genomePackage) || die "Could not create $genomePackage: $!";
            } else {
                $stats->Add(packageReplaced => 1);
                print "Clearing $genomePackage directory.\n";
                File::Copy::Recursive::pathempty($genomePackage);
            }
            # Create the metadata file.
            open(my $mh, ">$genomePackage/data.tbl") || die "Could not create metadata file in $genomePackage: $!";
            for my $key (sort keys %data) {
                print $mh "$key\t$data{$key}\n";
            }
            close $mh;
            # Copy the fasta file.
            File::Copy::Recursive::fcopy($contigFile, "$genomePackage/bin.fa") || die "Could not copy FASTA file into $genomePackage: $!";
            # Copy the GTO file.
            File::Copy::Recursive::fcopy($gtoFile, "$genomePackage/bin.gto") || die "Could not copy GTO file into $genomePackage: $!";
        }
    }
    $stats->Add(binProcessed => 1);
    $binNum++;
}
print "All done.\n" . $stats->Show();
