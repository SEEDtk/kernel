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


package Bin::Package;

    use strict;
    use warnings;
    use Bin;
    use GenomeTypeObject;
    use File::Copy::Recursive;

=head1 Genome Package Utilities

This module contains methods for manipulating genome packages.Each genome package is a directory with the
same name as the bin's genome ID and contains the following files

=over 4

=item bin.gto

A L<GenomeTypeObject> for the genome formed from the bin.

=item bin.fa

A FASTA file containing the bin's contigs.

=item data.tbl

A tab-delimited file of key/value pairs, with the following keys.

=over 8

=item Genome Name

Name of the bin genome

=item Sample Name

Name of the sample from which the bin came.

=item Bin Number

The ID number of the bin in the sample.

=item Contigs

The number of contigs in the bin.

=item Base pairs

The number of DNA base pairs in the bin's contigs.

=item Ref Genome

The ID of the closest reference genome.

=item Ref Name

The name of the closest reference genome.

=back

=back

For packages that do not come from bins, the C<Sample Name>, C<Bin Number>, C<Ref Genome>, and C<Ref Name>
fields will be missing, and other fields describing the package source will be substituted.

=head2 Public Methods

=head3 CreateFromSample

    Bin::Package::CreateFromSample($dir, $sampleName, $stats, $force, $outDir);

Create genome packages from a processed sample directory.

=over 4

=item dir

The full path of the sample directory.

=item sampleName

The sample name.

=item stats

A L<Stats> object for tracking activity.

=item force

TRUE if an existing package directory should be overwritten, else FALSE.

=item outDir

The name of the output directory containing the genome packages.

=back

=cut

sub CreateFromSample {
    my ($dir, $sampleName, $stats, $force, $outDir) = @_;
    # Read in the scores for the reference genomes. We use these to compute the closest genome to each bin.
    # The hash will map each genome ID to a [score, name] tuple.
    my %refGenomes;
    print "Processing reference genomes.\n";
    open(my $ih, '<', "$dir/ref.genomes.scores.tbl") || die "Could not open ref genomes score table: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /\t(\d+\.\d+)\t([\d\.]+)\t(.+)$/) {
            $refGenomes{$1} = [$2, $3];
        }
    }
    # Now read in the bins.
    print "Reading bin file.\n";
    my $binList = Bin::ReadBins("$dir/bins.rast.json");
    # Create a map from contig IDs to bins.
    my %binMap;
    for my $bin (@$binList) {
        $binMap{$bin->contig1} = $bin;
    }
    # Now we process the bins sequentially. For each GTO/FA pair, we search the contigs to find the matching bin object.
    # The bin object produces the majority of the data.tbl values.
    my $binN = 1;
    while (-f "$dir/bin$binN.gto") {
        my $binName = "bin$binN";
        print "Processing $binName.\n";
        $stats->Add(bins => 1);
        my $gto = GenomeTypeObject->create_from_file("$dir/$binName.gto");
        # Get the ID and name of the genome.
        my $genomeID = $gto->{id};
        my $name = $gto->{scientific_name};
        # Compute the output directory name.
        my $genomeDir = "$outDir/$genomeID";
        if (-d $genomeDir && ! $force) {
            print "Package already exists for $binName: $genomeID.\n";
            $stats->Add(binsAlreadyFound => 1);
        } else {
            # Here the bin is new or we are forcing.
            if (! -d $genomeDir) {
                print "Creating $genomeDir.\n";
                mkdir $genomeDir;
                $stats->Add(binProcessed => 1);
            } else {
                print "Replacing $genomeDir.\n";
                File::Copy::Recursive::pathempty($genomeDir);
                $stats->Add(binReplaced => 1);
            }
            # Find the bin object. One of the contigs will identify the bin. We stop when we hit it.
            my $bin;
            my $contigs = $gto->{contigs};
            for my $contig (@$contigs) { last if $bin;
                $bin = $binMap{$contig->{id}};
            }
            die "Bin object not found for $binName." if ! $bin;
            # Copy the main files.
            File::Copy::Recursive::fcopy("$dir/$binName.gto", "$genomeDir/bin.gto") ||
                die "Error copying $binName.gto: $!";
            File::Copy::Recursive::fcopy("$dir/$binName.fa", "$genomeDir/bin.fa") ||
                die "Error copying $binName.fa: $!";
            # We will compute the data table values in here.
            my %data;
            $data{'Genome Name'} = $name;
            $data{'Sample Name'} = $sampleName;
            $data{'Bin Number'} = $binN;
            $data{'Contigs'} = $bin->contigCount;
            $data{'Base pairs'} = $bin->len;
            # Find the closest reference genome.
            my @refs = $bin->refGenomes;
            my ($closest, $best) = ('', 0);
            for my $ref (@refs) {
                my $score = $refGenomes{$ref}[0];
                if ($score >= $best) {
                    $best = $score;
                    $closest = $ref;
                }
            }
            $data{'Ref Genome'} = $closest;
            $data{'Ref Name'} = $refGenomes{$closest}[1];
            # Write the data file.
            open(my $oh, '>', "$genomeDir/data.tbl") || die "Could not open $binName data table file: $!";
            for my $key (sort keys %data) {
                print $oh "$key\t$data{$key}\n";
            }
        }
        # Move to the next bin.
        $binN++;
    }
}

1;