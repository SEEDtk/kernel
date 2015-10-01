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


package Bin::Compute;

    use strict;
    use warnings;
    use Bin;
    use Stats;
    use IO::File;

=head1 Compute the Bins For a Community's Contigs.

This object contains the main method for converting contigs into bins. It takes as input a L<Bin::Score>
object containing scoring parameters, plus a list of bins, and an optional file containing scoring
vectors.

The fields of the object are as follows.

=over 4

=item logh

Open file handle for producing log and trace output.

=item stats

Statistics object for keeping statistics about a run.

=item score

L<Bin::Score> object for computing comparison scores between bins.

=back

=head2 Special Methods

    my $computer = Bin::Compute->new($score, %options);

Create a new bin computation object.

=over 4

=item score

A L<Bin::Score> object for comparing two bins.

=item options

A hash of optional parameters.

=over 8

=item logFile

An open file handle for the log/trace file, or the name of the log/trace file. The default is to write to STDOUT.

=back

=back

=cut

sub new {
    my ($class, $score, %options) = @_;
    # Handle the log file.
    my $logh;
    if (! $options{logFile}) {
        $logh = \*STDOUT;
    } elsif (! ref $options{logFile}) {
        $logh = IO::File->new(">$options{logFile}") || die "Could not open log: $!";
        $logh->autoflush(1);
    } else {
        $logh = $options{logFile};
    }
    # Create the object.
    my $retVal = {
        score => $score,
        logh => $logh,
        stats => Stats->new()
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Query Methods

=head3 stats

    my $stats = $helper->stats;

Return the L<Stats> object for tracking statistics.

=cut

sub stats {
    return $_[0]->{stats};
}


=head2 Public Methods

=head3 ProcessScores

    my $bins = $computer->ProcessScores($binList, $vectorFile);

Score comparisons between the specified bins to merge them into larger bins. We will first perform a comparison between
all the bins with found reference genomes. The non-zero comparisons will be sorted and used to combine the bins.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item vectorFile (optional)

If specified, a tab-delimited file containing comparison vectors for all the interesting contigs. For each pair of
interesting contigs, there will be a record containing (0) the first contig ID, (1) the second contig ID, (2) the
coverage score, (3) the tetranucleotide score, (4) the closest-reference-genome score, (4) the number of universal
roles not in common, and (5) the number of universal roles in common. If the file is not provided, the scores will
be computed from scratch.

=item RETURN

Returns a list of L<Bin> objects, sorted from most interesting to least, representing merged copies of the original
bins.

=back

=cut

sub ProcessScores {
    my ($self, $binList, $vectorFile) = @_;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Get the scoring object and the log file handle.
    my $score = $self->{score};
    my $logh = $self->{logh};
    # Get the number of bins.
    my $binCount = scalar @$binList;
    print $logh "$binCount input contigs specified for score processing.\n";
    print $logh $score->Show();
    $stats->Add(contigsRead => $binCount);
    my $filteredCount = scalar @$binList;
    my $totalScores = $filteredCount * ($filteredCount + 1) / 2;
    print $logh "$filteredCount contigs found in input. $totalScores scores required.\n";
    # This list contains 3-tuples for each pair of contigs containing the contig IDs and the comparison vector.
    my @scores;
    # This will count our progress.
    my $scoresComputed = 0;
    # Do we have a vector file?
    if ($vectorFile) {
        # Yes. Read from the file.
        print $logh "Reading score vectors from $vectorFile.\n";
        open(my $vh, "<", $vectorFile) || die "Could not open score vector file: $!";
        while (! eof $vh) {
            my $line = <$vh>;
            chomp $vh;
            my ($contigI, $contigJ, @vector) = split /\t/, $line;
            my $scoreValue = $score->ScoreV(\@vector);
            $stats->Add(pairsScored => 1);
            $scoresComputed++;
            if ($scoresComputed % 100000 == 0) {
                print $logh "$scoresComputed of $totalScores scores computed.\n";
            }
            if ($scoreValue > 0) {
                $stats->Add(pairsKept => 1);
                push @scores, [$contigI, $contigJ, $scoreValue];
            }
        }
    } else {
        # No vector file. Compute the scores manually.
        my ($i, $j);
        # Now the nested loop.
        for ($i = 0; $i < $filteredCount; $i++) {
            my $binI = $binList->[$i];
            my $contigI = $binI->contig1;
            for ($j = $i + 1; $j < $filteredCount; $j++) {
                my $binJ = $binList->[$j];
                my $contigJ = $binJ->contig1;
                my $scoreValue = $score->Score($binI, $binJ);
                $scoresComputed++;
                if ($scoresComputed % 1000000 == 0) {
                    print $logh "$scoresComputed of $totalScores scores computed.\n";
                }
                $stats->Add(pairsScored => 1);
                if ($scoreValue > 0) {
                    $stats->Add(pairsKept => 1);
                    push @scores, [$contigI, $contigJ, $scoreValue];
                }
            }
        }
    }
    # Use the score vector to cluster the contigs.
    my $retVal = $self->ClusterContigs($binList, \@scores);
    return $retVal;
}


=head3 ClusterContigs

    my $bins = $computer->ClusterContigs($binList, $scores, %options);

Cluster contigs into bins based on a score list. Essentially, any two contigs that should go together
should have a score in the score list. We attempt to make the largest bins we can out of the scored
connections.

Two multi-contig bins are normally re-scored to see if they should be combined. If the C<closed> option
is specified, there is no rescoring and bins are always combined.

=over 4

=item binList

A reference to a list of L<Bin> objects containing the input contigs, one per bin.

=item scores

A reference to a list of 3-tuples, each tuple consisting of (0) a first contig ID, (1) a second contig ID for
a contig that may belong with the first, and (2) a score for the combination.

=item options

A hash of options.

=over 8

=item closed

If TRUE, then bins are not re-scored before being combined. The incoming pairs are assumed to imply transitivity.

=back

=item RETURN

Returns a list of bins in order from most to least interesting.

=back

=cut

sub ClusterContigs {
    my ($self, $binList, $scores, %options) = @_;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Get the scoring object and the log file handle.
    my $score = $self->{score};
    my $logh = $self->{logh};
    # Get the option on whether or not we re-score.
    my $open = ! $options{closed};
    # We must loop through the contigs, comparing them. This hash tells us which bin
    # contains each contig.
    my %contig2Bin = map { $_->contig1 => $_ } @$binList;
    print $logh "Sorting " . scalar(@$scores) . " scores.\n";
    my @scores = sort { $b->[2] <=> $a->[2] } @$scores;
    print $logh "Merging bins.\n";
    for my $scoreTuple (@scores) {
        # Get the bins.
        my ($contigI, $contigJ, $scoreValue) = @$scoreTuple;
        my $binI = $contig2Bin{$contigI};
        my $binJ = $contig2Bin{$contigJ};
        # Compute the score for these two bins.
        if (! $binI) {
            print $logh "WARNING: Missing bin for $contigI.\n";
            $scoreValue = 0;
            $stats->Add(missingBin => 1);
        } elsif (! $binJ) {
            print $logh "WARNING: Missing bin for $contigJ.\n";
            $scoreValue = 0;
            $stats->Add(missingBin => 1);
        } elsif ($binI->contig1 eq $binJ->contig1) {
            # Here the contigs are already in the same bin.
            $scoreValue = 0;
            $stats->Add(binsAlreadyMerged => 1);
        } elsif ($open && ($binI->contigCount > 1 || $binJ->contigCount > 1)) {
            # The bins are not singletons. We have to re-compute the score.
            $stats->Add(complexBinMatch => 1);
            $scoreValue = $score->Score($binI, $binJ);
        }
        # If the score is nonzero, we can merge them.
        if ($scoreValue > 0) {
            $stats->Add(binsMerged => 1);
            # The bins will be merged into bin I. First, we must update the contig-to-bin map.
            for my $contigJX ($binJ->contigs) {
                $contig2Bin{$contigJX} = $binI;
            }
            # Now merge the bins.
            $binI->Merge($binJ);
        }
    }
    # Get the final list of bins.
    print $logh "Preparing for output.\n";
    my %found;
    my @bins;
    for my $contigX (keys %contig2Bin) {
        my $bin = $contig2Bin{$contigX};
        my $binID = $bin->contig1;
        if (! $found{$binID}) {
            push @bins, $bin;
            $found{$binID} = 1;
            $stats->Add(outputBin => 1);
        }
    }
    print $logh scalar(@bins) . " bins output.\n";
    my @sortedBins = sort { Bin::cmp($a, $b) } @bins;
    # Return the statistics and the list of bins.
    return \@sortedBins;
}

=head3 ProcessBlast

    my $stats = Bin::Compute::ProcessBlast($shrub, $workDir, \@refGenomes, \%contigBins, $contigFasta, \@roles, %options);

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
                # If this is a universal role, count it.
                if ($functionID) {
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