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


package Bin::Blast;

    use strict;
    use warnings;
    use gjo::BlastInterface;
    use SeedUtils;
    use Stats;
    use File::Spec;
    use File::Copy::Recursive;

=head1 Blast Analysis Object

This library creates a BLAST database for a community sample, then uses it to find the closest reference genomes
to each contig in the sample. It accepts as input a selected list of reference genomes to use. It is presumed
a prior process has been used to select reference genomes likely to produce a hit. The method is
also responsible for returning the universal roles found in each contig.

Because there are so many more contigs than reference genomes, a BLAST database is created for the contigs and
then each reference genome is blasted against it. Because the contigs are DNA and we are blasting proteins, the
tool used will be C<tblastn>.

When the query sequences are created, the comment field for each sequence will be the ID of its functional
assignment. This can be used to determine whether or not a universal role has been hit.

=head2 Public Methods

=head3 Process

    my $stats = Bin::Blast::Process($shrub, $workDir, \@refGenomes, \%contigBins, $contigFasta, \@roles, %options);

BLAST a set of contigs and update the L<Bin> objects for those contigs.

=over 4

=item shrub

The L<Shrub> object for accessing the database.

=item workDir

Name of the working directory to contain intermediate files.

=item refGenomes

Reference to a list of reference genome IDs.

=item contigBins

Reference to a hash mapping each sample contig ID to a L<Bin> object which is to contain the results of the blasting.

=item contigFasta

Name of the FASTA file containing the sample contigs.

=item roles

Reference to a list of universal role IDs.

=item options

A hash of options. These include

=over 8

=item minsim

The minimum acceptable percent identity similarity. The default is C<0.9>.

=item minlen

The minimum acceptable match length. The default is C<150>.

=item priv

Privilege level for the functions of the feature. The default is C<1>.

=item maxE

The maximum accetable E-value. The default is C<1e-50>.

=back

=item RETURN

Returns a statistics object containing useful information about the BLAST results.

=back

=cut

sub Process {
    my ($shrub, $workDir, $refGenomes, $contigBins, $contigFasta, $roles, %options) = @_;
    # This will contain the return statistics.
    my $stats = Stats->new();
    # Get the options.
    my $minsim = $options{minsim} // 0.9;
    my $minlen = $options{minlen} // 150;
    my $priv = $options{priv} // 1;
    my $maxE = $options{maxE} // 1e-50;
    # Create a hash of the universal role IDs.
    my %uniRoles = map { $_ => 1 } @$roles;
    # Insure we have a copy of the sample contigs in the working directory. Note we need to deal with
    # trailing-slash craziness.
    my $blastFasta = $contigFasta;
    my $absWorkDir = File::Spec->rel2abs($workDir);
    my $absContigFasta = File::Spec->rel2abs($contigFasta);
    my ($volWorkDir, $pathWorkDir) = File::Spec->splitpath($absWorkDir, 1);
    my ($volFasta, $pathFasta, $nameFasta) = File::Spec->splitpath($absContigFasta);
    $pathFasta =~ s/[\/\\]$//;
    $pathWorkDir =~ s/[\/\\]$//;
    if ($pathFasta ne $pathWorkDir) {
        my $newName = "$workDir/$nameFasta";
        print "Copying $contigFasta to $newName.\n";
        File::Copy::Recursive($contigFasta, $newName);
        $blastFasta = $newName;
    }
    # Create the BLAST database for the sample contigs. If it already exists, it will be reused.
    my $blastDbName = gjo::BlastInterface::get_db($blastFasta, 'tblastn', $workDir);
    print "BLAST database found at $blastDbName.\n";
    # The contig bins will contain the universal role information. We cannot, however, track the closest
    # genomes there because we want only the best two. This hash will map each sample contig ID to a
    # sub-hash of genome scores. For each genome the sub-hash will contain the percent identity score of its best hit.
    my %contigGenomes;
    my $totalGenomes = scalar (@$refGenomes);
    # Loop through the reference genomes.
    for my $refGenome (@$refGenomes) {
        my $blasted = $stats->Add(refGenomesBlasted => 1);
        print "Processing $refGenome for BLAST ($blasted of $totalGenomes).\n";
        # Create the FASTA file for this genome's query sequence.
        my $queryFileName = "$workDir/$refGenome.fa";
        open(my $oh, ">", $queryFileName) || die "Could not create FASTA output file: $!";
        # Get the sequences and functions for all of this genome's proteins.
        my @tuples = $shrub->GetAll('Feature Protein AND Feature Feature2Function', 'Feature(id) LIKE ? AND Feature2Function(security) = ?',
            ["fig|$refGenome.peg.%", $priv], 'Feature(id) Feature2Function(to-link) Protein(sequence)');
        # Write them in FASTA format.
        for my $tuple (@tuples) {
            my ($id, $function, $sequence) = @$tuple;
            # Only keep the function if it is a universal role.
            if (! $uniRoles{$function}) {
                $function = '';
            } else {
                $stats->Add(uniRoleProteins => 1);
            }
            print $oh ">$id $function\n$sequence\n";
            $stats->Add(genomeProteins => 1);
        }
        close $oh;
        # Blast this genome against the sample contigs.
        my $matches = gjo::BlastInterface::blast($queryFileName, $blastDbName, 'tblastn',
            { outForm => 'hsp', minIden => $minsim, minLen => $minlen, maxE => $maxE });
        my $matchCount = scalar @$matches;
        $stats->Add(blastMatches => $matchCount);
        if ($matchCount) {
            print "$matchCount hits found.\n";
            # Loop through the matches. Note they are pre-filtered for length and percent identity.
            for my $match (@$matches) {
                # Get the pieces of the HSP object.
                my $functionID = $match->[1];
                my $contigID = $match->[3];
                my $score = $match->[11] / $match->[10];
                # Check to see if this is the genome's best score for this contig.
                my $oldScore = $contigGenomes{$contigID}{$refGenome} // 0;
                if ($oldScore < $score) {
                    $contigGenomes{$contigID}{$refGenome} = $score;
                }
                # If this is a universal role and we match over 60% of the length, count it.
                if ($functionID && $match->[10] >= 0.6 * $match->[2]) {
                    $stats->Add(uniRoleFound => 1);
                    $contigBins->{$contigID}->add_prots($functionID);
                }
            }
        }
    }
    # Now all the reference genomes have been blasted. For each contig, we need to choose the best two.
    # Note that many contigs will not have any hits. These do not even appear in the hash, so we don't
    # worry about skipping them.
    print "Assigning genomes to sample contigs.\n";
    for my $contigID (keys %contigGenomes) {
        $stats->Add(contigsWithRefGenomes => 1);
        my $genomeH = $contigGenomes{$contigID};
        my @sorted = sort { $genomeH->{$b} <=> $genomeH->{$a} } keys %$genomeH;
        my $maxN = scalar @sorted;
        my @genomesL;
        # This is a bit tricky. If multiple genomes have the same score, we keep them even if it means
        # we are passing back more than 2.
        if ($maxN <= 2) {
            @genomesL = @sorted;
        } else {
            @genomesL = shift @sorted;
            my $g2 = shift @sorted;
            push @genomesL, $g2;
            my $g2score = $genomeH->{$g2};
            while (scalar(@sorted) && $genomeH->{$sorted[0]} == $g2score) {
                push @genomesL, shift @sorted;
            }
        }
        # Add the genomes found to the contig's bin.
        $contigBins->{$contigID}->add_ref(@genomesL);
    }
    # Return the statistics.
    return $stats;
}


1;