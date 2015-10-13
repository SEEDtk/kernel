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
    use BasicLocation;

=head1 Blast Analysis Object

This object creates a BLAST database for a community sample. Several methods are provided for doing different
types of BLASTing against it. The query sequences will be proteins, and the sample is DNA, so the
tool used will be C<tblastn>.

When the query sequences are created, the comment field for each sequence will be the ID of its functional
assignment. This can be used to determine whether or not a universal role has been hit.

The fields of this object are as follows.

=over 4

=item uniRolesH

Reference to a hash that maps each universal role ID to its description.

=item shrub

Reference to a L<Shrub> object for accessing the database.

=item workDir

Name of the working directory to contain intermediate and data files.

=item blastDb

Name of a BLAST database containing the sample contigs.

=item priv

Privilege level for computing the functions of the feature. The default is C<1>.

=item maxE

The maximum acceptable E-value. The default is C<1e-5>.

=item defaultRole

The default universal role.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<100>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>).

=back

=head2 Special Methods

=head3 new

    my $blaster = Bin::Blast->new($shrub, $workDir, $contigFasta, %options);

Create a BLAST database for the sample contigs.

=over 4

=item shrub

The L<Shrub> object for accessing the database.

=item workDir

Name of the working directory to contain intermediate files.

=item contigFasta

Name of the FASTA file containing the sample contigs.

=item options

A hash of options. These include

=over 8

=item priv

Privilege level for computing the functions of the feature. The default is C<1>.

=item uniRoles

Name of a tab-delimited file containing the universal roles. Each record must contain (0) a universal role ID, (1) the
number of times it was found in the database, and (2) the role's description. The default is C<uni_roles.tbl> in the
global data directory.

=item maxE

The maximum acceptable E-value. The default is C<1e-50>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<100>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=back

=back

=cut

sub new {
    my ($class, $shrub, $workDir, $contigFasta, %options) = @_;
    # Get the options.
    my $priv = $options{priv} // 1;
    my $maxE = $options{maxE} // 1e-5;
    my $gap = $options{gap} // 100;
    my $minlen = $options{minlen} // 0.50;
    # This will track the best universal role.
    my ($defaultRole, $drCount) = ('', 0);
    # Create a hash of the universal roles.
    my %uniRoleH;
    my $uniRoleFile = $options{uniRoles} // "$FIG_Config::global/uni_roles.tbl";
    open(my $ih, "<$uniRoleFile") || die "Could not open universal role file: $!";
    while (! eof $ih) {
        my ($role, $count, $desc) = SeedUtils::fields_of($ih);
        $uniRoleH{$role} = $desc;
        if ($count > $drCount) {
            ($defaultRole, $drCount) = ($role, $count);
        }
    }
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
        File::Copy::Recursive::fcopy($contigFasta, $newName);
        $blastFasta = $newName;
    }
    # Create the BLAST database for the sample contigs. If it already exists, it will be reused.
    my $blastDbName = gjo::BlastInterface::get_db($blastFasta, 'tblastn', $workDir);
    # Create the object.
    my $retVal = {
        priv => $priv,
        maxE => $maxE,
        uniRoleH => \%uniRoleH,
        shrub => $shrub,
        blastDb => $blastDbName,
        workDir => $workDir,
        defaultRole => $defaultRole,
        gap => $gap,
        minlen => $minlen
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Query Methods

=head3 uni_hash

    my $uni_hash = $blaster->uni_hash;

Return the universal role hash.

=cut

sub uni_hash {
    my ($self) = @_;
    return $self->{uniRoleH};
}

=head3 uni_total

    my $uni_total = $blaster->uni_total;

Return the total number of universal roles.

=cut

sub uni_total {
    my ($self) = @_;
    return scalar keys %{$self->{uniRoleH}};
}


=head2 Public Manipulation Methods

=head3 BestProteinHits

    my $hitHash = $blaster->BestProteinHits($refGenome, $binFasta);

Find the best hit in the specified bin for each universal protein in the specified reference genome.

=over 4

=item refGenome

ID of the reference genome whose universal proteins are to be analyzed.

=item binFasta

Name of a DNA FASTA file for the bin's contigs.

=item RETURN

Returns a reference to a hash mapping each universal protein to its best hit in the bin.

=back

=cut

sub BestProteinHits {
    my ($self, $refGenome, $binFasta) = @_;
    # Declare the return variable.
    my %retVal;
    # Get the options.
    my $gap = $self->{gap};
    # Get the universal protein FASTA.
    my ($queryFileName, $uniLens) = $self->GetRefGenomeFasta($refGenome, uniOnly => 1);
    # Get the list of hits.
    my @matches = sort { ($a->sid cmp $b->sid) or ($a->s1 <=> $b->s1) }
            gjo::BlastInterface::blast($queryFileName, $binFasta, 'tblastn',
            { outForm => 'hsp', maxE => $self->{maxE} });
    # Now we track the matches for each universal protein. The matches come back in the form of Hsp objects.
    # We want the match with the best percentage (defined as the identity count over the match length). If two
    # matches have the same percentage, we pick the longer one. We also need to merge matches to compensate for
    # frame shifts. In this case, the percentages are combined. Our approach will be to do the merging and then
    # throw out all but the best for each protein. This hash lists the matches for each protein. Each match is
    # stored in here in the form [location, n_id];
    my %uniMatches;
    for my $match (@matches) {
        my $prot = $match->qdef;
        if (! exists $uniMatches{$prot}) {
            $uniMatches{$prot} = [[$match->sloc, $match->n_id]];
        } else {
            my $protMatches = $uniMatches{$prot};
            # Now we have a list of the matches for this protein. If the new match belongs with an old one, it will
            # belong with the last.
            my $oldMatch = pop @$protMatches;
            my ($oldLoc, $old_n_id) = @$oldMatch;
            my $newLoc = $match->sloc;
            if ($oldLoc->Dir eq $newLoc->Dir && ($newLoc->Left - $oldLoc->Right) < $self->{gap} ) {
                # We belong to the previous hit.
                $oldLoc->Merge($newLoc);
                $old_n_id += $match->n_id;
                push @$protMatches, [$oldLoc, $old_n_id];
            } else {
                # We are a new hit.
                push @$protMatches, $oldMatch, [$newLoc, $match->n_id];
            }
        }
    }
    # Now select the best hit for each protein and put it in the return hash. Only hits greater than the minimum length are
    # even considered.
    for my $prot (keys %uniMatches) {
        my $protMatches = $uniMatches{$prot};
        my ($bestMatch, $bestScore) = (undef, 0);
        for my $protMatch (@$protMatches) {
            my ($matchLoc, $match_n_id) = @$protMatch;
            # Compute this match's score.
            my $newScore = $match_n_id / $matchLoc->Length;
            if ($newScore > $bestScore || $newScore == $bestScore && $matchLoc->Length > $bestMatch->Length) {
                # Here we have a new best.
                $bestMatch = $matchLoc;
                $bestScore = $newScore;
            }
        }
        if ($bestMatch) {
            $retVal{$prot} = $bestMatch;
        }
    }
    # Return the hash of best matches.
    return \%retVal;
}

=head3 BestDnaHits

    my $bestHitHash = $blaster->BestDnaHits(\%seqHash, $fastaFile);

Find the best protein hits for a set of DNA sequences. The DNA sequences will be passed in via a hash that maps IDs to
sequences. The output hash will return, for each ID, the ID and definition string for the protein in the specified FASTA
file that represents the corresponding DNA sequence's best hit.

=over 4

=item seqHash

Reference to a hash that maps identifiers to DNA sequences.

=item fastaFile

The name of a FASTA file containing protein sequences.

=item RETURN

Returns a hash mapping each sequence identifier to a 2-tuple consisting of (0) the id and (1) the definition string of the
best-hit protein.

=back

=cut

sub BestDnaHits {
    my ($self, $seqHash, $fastaFile) = @_;
    # Create a FASTA triples object for use by the blast interface software.
    my @triples = map { [$_, '', $seqHash->{$_}] } keys %$seqHash;
    # BLAST to get the matches.
    my $matches = gjo::BlastInterface::blast(\@triples, $fastaFile, 'blastx',
            { maxE => $self->{maxE}, tmp_dir => $self->{workDir}, outForm => 'hsp' });
    # The matches come back in the order of best match to worst. For each match, the query sequence
    # ID will be the identifier from $seqHash, and the subject sequence ID will be a protein ID. The
    # first hit found for each query sequence will be the best.
    my %retVal;
    for my $match (@$matches) {
        my $seqID = $match->qid;
        if (! exists $retVal{$seqID}) {
            my $found = [$match->sid, $match->sdef // ''];
            $retVal{$seqID} = $found;
        }
    }
    # Return the result hash.
    return \%retVal;
}


=head3 FilterBestHits

    my $bbhitHash = $blaster->FilterBestHits(\%hitHash, $refGenome);

This method filters out the hits from L<BestProteinHits> and only keeps the bidirectional best hits.
Each protein is associated with a location in the bin. We blast that DNA's location back against the complete
set of proteins in the genome and verify that the first hit for each location (the best) is the correct protein.

=over 4

=item hitHash

Reference to a hash mapping protein role IDs to their best hit sequences.

=item refGenome

ID of the reference genome against which the sequences are to be blasted.

=item RETURN

Returns a reference to a hash whose keys are the universal roles whose best hits are bidirectional.

=back

=cut

sub FilterBestHits {
    my ($self, $hitHash, $refGenome) = @_;
    # Get the working directory.
    my $workDir = $self->{workDir};
    # Get the shrub database object.
    my $shrub = $self->{shrub};
    # Create a FASTA file for the reference genome.
    my ($fastaFileName) = $self->GetRefGenomeFasta($refGenome);
    # This will be the return hash.
    my %retVal;
    # Perform the DNA blast.
    my $bestDnaHits = $self->BestDnaHits($hitHash, $fastaFileName);
    # Return the DNA hits back to the same protein.
    for my $prot (keys %$bestDnaHits) {
        my ($hitID, $hitProt) = @{$bestDnaHits->{$prot}};
        if ($hitProt eq $prot) {
            $retVal{$prot} = 1;
        }
    }
    # Return the result.
    return \%retVal;
}

=head3 FindProtein

    my $hitHash = $blaster->FindProtein(\@sampleGenomes, $funID);

Return a list of hits representing the best occurrences of a particular role in the sample contigs.
At most one hit per contig will be returned.

=over 4

=item sampleGenomes

Reference to a list of the IDs of the genomes from which to take the prototype versions of the role. If
omitted, C<83333.1> will be assumed.

=item funID

ID of the desired function. If omitted, the default (most common) universal role will be assumed.

=item RETURN

Returns a reference to a hash mapping each contig ID for which there was a hit to a L<BasicLocation> object
describing the location that is the best match for the specified protein.

=back

=cut

sub FindProtein {
    my ($self, $sampleGenomes, $funID, %options) = @_;
    # Declare the return variable.
    my %retVal;
    # Get the options.
    my $minlen = $self->{minlen};
    # Process the defaults.
    $sampleGenomes //= ['83333.1'];
    $funID //= $self->{defaultRole};
    # Insure we have a list reference for the genomes.
    if (! ref $sampleGenomes) {
        $sampleGenomes = [$sampleGenomes];
    }
    # Get the query protein sequences. Note that the query will return a set of FASTA triples such as would be
    # expected by BlastInterface.
    my @in = map { '?' } @$sampleGenomes;
    my $filter = "(" . join(", ", @in) . ")";
    my (@prots) = $self->{shrub}->GetAll('Function2Feature Feature Feature2Genome AND Feature Protein',
            "Function2Feature(from-link) = ? AND Function2Feature(security) = ? AND Feature2Genome(to-link) IN $filter",
            [$funID, $self->{priv}, @$sampleGenomes], 'Feature(id) Function2Feature(from-link) Protein(sequence)');
    if (! @prots) {
        die "Protein $funID not found in any sample genomes.";
    }
    # Compute the minimum match length by taking the mean of all the proteins found.
    my ($protLen, $protCount) = (0, 0);
    for my $prot (@prots) {
        my $len = length($prot->[2]);
        $protLen += $len;
        $protCount++;
    }
    my $minDnaLen = int($minlen * $protLen / $protCount) * 3;
    # Look for matches.
    my @matches = sort { ($a->sid cmp $b->sid) or ($a->s1 <=> $b->s1) }
            gjo::BlastInterface::blast(\@prots, $self->{blastDb}, 'tblastn',
            { outForm => 'hsp', maxE => $self->{maxE} });
    # The matches are in the form of Hsp objects. They are sorted by start position within contig.
    # We condense the matches into location objects in the following hash. This is where the gap
    # resolution is processed.
    my %contigs;
    for my $match (@matches) {
        $self->MergeMatch(\%contigs, $match);
    }
    # For each contig, find the largest hit location that is greater than the computed minimum DNA length.
    # This will be put in the return hash.
    for my $contig (keys %contigs) {
        my ($match) = sort { $b->Length <=> $a->Length} @{$contigs{$contig}};
        if ($match->Length >= $minDnaLen) {
            $retVal{$contig} = $match;
        }
    }
    # Return the result.
    return \%retVal;
}

=head3 MatchProteins

    my $contigHash = $blaster->MatchProteins($seqHash, $funID, $count);

BLAST contig sequences against all occurrences of a protein in the Shrub database. The protein is identified by its
functional assignment ID.

=over 4

=item seqHash

Reference to a hash mapping each contig ID to a DNA sequence.

=item funID

ID of the desired function. If omitted, the default universal role is used.

=item count

Number of reference genomes to keep for each bin. If omitted, the default is C<1>.

=item RETURN

Returns a reference to a hash mapping each incoming contig ID to a reference to a list of reference genome IDs.

=back

=cut

sub MatchProteins {
    my ($self, $seqHash, $funID, $count) = @_;
    # Get the defaults.
    $funID //= $self->{defaultRole};
    $count //= 1;
    # Get the database.
    my $shrub = $self->{shrub};
    # Create the BLAST database for the protein.
    my $protFastaFile = $self->{workDir} . "/$funID.fa";
    if (! -s $protFastaFile) {
        # Here it does not exist, so we must create it.
        $self->ProtDnaDatabase($funID, $protFastaFile);
    }
    # Create a list of FASTA triplets from the matched sequences.
    my @triples = map { [$_, '', $seqHash->{$_}] } keys %$seqHash;
    # BLAST to get the matches against the protein DNA from the database.
    my $matches = gjo::BlastInterface::blast(\@triples, $protFastaFile, 'blastn',
            { maxE => $self->{maxE}, tmp_dir => $self->{workDir}, outForm => 'hsp' });
    # The query sequence IDs are sample contigs representing bins. The subject sequence IDs are
    # genome IDs. Save the best N genomes for each sample contig.
    my %retVal;
    for my $match (@$matches) {
        my ($contig, $genome) = ($match->qid, $match->sid);
        if (! $retVal{$contig}) {
            $retVal{$contig} = [$genome];
        } else {
            my $gList = $retVal{$contig};
            # If this is a new genome and we do not already have too many, keep it.
            if (! (grep { $_ eq $genome } @$gList) && scalar(@$gList) < $count) {
                push @$gList, $genome;
            }
        }
    }
    # Return the hash of contig IDs to genome information.
    return \%retVal;
}

=head3 Process

    my $stats = $blaster->Process(\%contigBins, \@refGenomes);

BLAST the specified reference genomes against the contigs and assign the best reference genome to each contig.

=over 4

=item contigBins

Reference to a hash mapping each sample contig ID to a L<Bin> object which is to contain the results of the blasting.

=item refGenomes

Reference to a list of reference genome IDs.

=item RETURN

Returns a statistics object describing the results of the BLASTing.

=back

=cut

sub Process {
    my ($self, $contigBins, $refGenomes) = @_;
    # This will contain the return statistics.
    my $stats = Stats->new();
    # Get the working directory.
    my $workDir = $self->{workDir};
    # Get the universal role hash.
    my $uniRoles = $self->{uniRoleH};
    # Get the shrub database object.
    my $shrub = $self->{shrub};
    # Get the options.
    my $maxE = $self->{maxE};
    my $priv = $self->{priv};
    # The contig bins will contain the universal role information. We cannot, however, track the closest
    # genomes there because we want only the best one. This hash will map each sample contig ID to a
    # 2-tuple of [genome ID, score].
    my %contigGenomes;
    my $totalGenomes = scalar (@$refGenomes);
    # Loop through the reference genomes.
    for my $refGenome (@$refGenomes) {
        my $blasted = $stats->Add(refGenomesBlasted => 1);
        print "Processing $refGenome for BLAST ($blasted of $totalGenomes).\n";
        # This hash will track the minimum match length of each universal role.
        my ($queryFileName, $uniLens) = $self->GetRefGenomeFasta($refGenome);
        print scalar(keys %$uniLens) . " universal roles found in $refGenome.\n";
        # Blast this genome against the sample contigs.
        my $matches = gjo::BlastInterface::blast($queryFileName, $self->{blastDb}, 'tblastn',
            { outForm => 'hsp', maxE => $maxE });
        my $matchCount = scalar @$matches;
        $stats->Add(blastMatches => $matchCount);
        if ($matchCount) {
            print "$matchCount hits found.\n";
            # This hash will track universal role hits. It is a double hash keyed on universal role ID followed by
            # contig ID and maps to location objects.
            my %uniHits = map { $_ => {} } keys %$uniLens;
            # Loop through the matches.
            for my $match (sort { ($a->sid cmp $b->sid) or ($a->s1 <=> $b->s1) } @$matches) {
                # Get the pieces of the HSP object.
                my $functionID = $match->qdef;
                my $contigID = $match->sid;
                # The score is percent identity first, then length of match.
                my $score = [$match->n_id / $match->n_mat, $match->n_mat];
                # Check to see if this is the genome's best score for this contig.
                my $oldScore = [0, 0];
                if ($contigGenomes{$contigID}) {
                    $oldScore = $contigGenomes{$contigID}[1];
                }
                if ($oldScore->[0] < $score->[0] || $oldScore->[0] == $score->[0] && $oldScore->[1] < $score->[1]) {
                    $contigGenomes{$contigID} = [$refGenome, $score];
                }
                # If this is a universal role, merge it into the hash.
                if ($functionID) {
                    my $uniSubHash = $uniHits{$functionID};
                    $stats->Add(uniRoleFound => 1);
                    if ($self->MergeMatch($uniSubHash, $match)) {
                        $stats->Add(uniRoleMerged => 1);
                    }
                }
            }
            # Check for any universal role matches of sufficient length.
            my $roleMatches = 0;
            for my $role (keys %uniHits) {
                my $contigHits = $uniHits{$role};
                my $minLen = $uniLens->{$role};
                for my $contig (keys %$contigHits) {
                    my $hitList = $contigHits->{$contig};
                    for my $hit (@$hitList) {
                        if ($hit->Length >= $minLen) {
                            $contigBins->{$contig}->merge_prots($role);
                            $stats->Add(uniRoleAssigned => 1);
                            $roleMatches++;
                        }
                    }
                }
            }
            print "$roleMatches universal role hits found by $refGenome BLAST.\n";
        }
    }
    # Now all the reference genomes have been blasted. For each contig, we need to assign the reference
    # genome and any universal roles.
    print "Assigning genomes to sample contigs.\n";
    for my $contigID (keys %contigGenomes) {
        $stats->Add(contigsWithRefGenomes => 1);
        my $genome = $contigGenomes{$contigID}[0];
        # Add the genome found to the contig's bin.
        $contigBins->{$contigID}->add_ref($genome);
    }
    # Return the statistics.
    return $stats;
}


=head2 Internal Utility Methods

=head3 ProtDnaDatabase

    $blaster->ProtDnaDatabase($funID, $fileName);

Create a FASTA file of all the DNA sequences for the specified function. The function is presumed to be a universal protein, so
we expect only one occurrence in each genome having a single location. The code is optimized for this assumption.

=over 4

=item funID

ID of the function for the protein of interest.

=item fileName

Name of the FASTA output file.

=back

=cut

sub ProtDnaDatabase {
    my ($self, $funID, $fileName) = @_;
    # Get the shrub database.
    my $shrub = $self->{shrub};
    # Get the DNA repo directory.
    my $dnaDir = $shrub->DNArepo();
    # Get the function privilege level.
    my $priv = $self->{priv};
    # Get all the occurrences of the function in the database.
    my @funTuples = $shrub->GetAll('Function2Feature Feature Feature2Contig AND Feature Genome',
            'Function2Feature(from-link) = ? AND Function2Feature(security) = ? ORDER BY Feature2Contig(from-link), Feature2Contig(ordinal)',
            [$funID, $priv], 'Feature2Contig(from-link) Genome(id) Genome(contig-file) Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(dir) Feature2Contig(len)');
    # This hash will contain the DNA sequences found for each genome.
    my %genomeDNA;
    # This hash contains the feature ID for each genome. It is used to detect multi-occurring proteins.
    my %genomeF;
    # Loop through the tuples.
    for my $tuple (@funTuples) {
        my ($fid, $genome, $dnaFile, $contig, $beg, $dir, $len) = @$tuple;
        # Is this feature a duplicate?
        if (! $genomeF{$genome} || $genomeF{$genome} eq $fid) {
            # No. Save the feature ID.
            $genomeF{$genome} = $fid;
            # Open the DNA file.
            open(my $fh, "<$dnaDir/$dnaFile") || die "Could not open FASTA file for $genome: $!";
            # Find the contig.
            my $found;
            while (! eof $fh && ! $found) {
                my $line = <$fh>;
                if ($line =~ />(\S+)/) {
                    $found = ($1 eq $contig);
                }
            }
            # At this point, we should have just read the contig header.
            if (! $found) {
                die "Could not find $contig in $genome FASTA file $dnaFile.";
            }
            my @dna;
            while (! eof $fh && $found) {
                my $line = <$fh>;
                chomp $line;
                if (substr($line, 0, 1) eq '>') {
                    $found = 0;
                } else {
                    push @dna, $line;
                }
            }
            # Get the DNA and adjust for direction.
            my $dna = substr(join("", @dna), $beg-1, $len);
            if ($dir eq '-') {
                SeedUtils::rev_comp(\$dna);
            }
            # Save it.
            if (exists $genomeDNA{$genome}) {
                $genomeDNA{$genome} .= $dna;
            } else {
                $genomeDNA{$genome} = $dna;
            }
        }
    }
    # Create the output file.
    open(my $ofh, ">$fileName") || die "Could not open FASTA output file: $!";
    for my $genome (sort keys %genomeDNA) {
        print $ofh ">$genome\n";
        print $ofh "$genomeDNA{$genome}\n";
    }
}

=head3 GetRefGenomeFasta

    my ($queryFileName, \%uniLens) = $blaster->GetRefGenomeFasta($refGenome, %options);

Create a FASTA file containing all the proteins in a reference genome. The universal roles will have comments
identifying which universal role is assigned to a protein.

=over 4

=item refGenome

ID of the target reference genome.

=item options

A hash containing zero or more of the following keys.

=over 8

=item uniOnly

If TRUE, only proteins with universal roles will be included.

=back

=item RETURN

Returns a two-element list consisting of (0) the name of the file created and (1) a reference to a hash mapping each
universal role ID to the minimum length required to match it.

=back

=cut

sub GetRefGenomeFasta {
    my ($self, $refGenome, %options) = @_;
    # Get the working directory.
    my $workDir = $self->{workDir};
    # Get the shrub database object.
    my $shrub = $self->{shrub};
    # Get the universal role hash.
    my $uniRoles = $self->{uniRoleH};
    # Get the options.
    my $uniOnly = $options{uniOnly} // 0;
    # Create the FASTA file for this genome's query sequence.
    # Get the options.
    my $priv = $self->{priv};
    my $uniMod = ($uniOnly ? '.uni' : '');
    my $queryFileName = "$workDir/$refGenome$uniMod.fa";
    open(my $oh, ">", $queryFileName) || die "Could not create FASTA output file: $!";
    # This hash will track the minimum match length for each universal role.
    my %uniLens;
    # Get the sequences and functions for all of this genome's proteins.
    my @tuples = $shrub->GetAll('Feature Protein AND Feature Feature2Function', 'Feature(id) LIKE ? AND Feature2Function(security) = ?',
        ["fig|$refGenome.peg.%", $priv], 'Feature(id) Feature2Function(to-link) Protein(sequence)');
    # Write them in FASTA format.
    for my $tuple (@tuples) {
        my ($id, $function, $sequence) = @$tuple;
        # Only keep the function if it is a universal role.
        if (! $uniRoles->{$function}) {
            $function = '';
        } else {
            # Track the universal role match length.
            $uniLens{$function} = int(length($sequence) * $self->{minlen}) * 3;
        }
        if ($function || ! $uniOnly) {
            print $oh ">$id $function\n$sequence\n";
        }
    }
    close $oh;
    return ($queryFileName, \%uniLens);
}

=head3 MergeMatch

    my $flag = $blaster->MergeMatch(\%contigs, $match);

Store the location of a hit against a contig. The incoming hash should contain a list of hit locations for a specific
protein, keyed by contig ID. The incoming match should be located to the right of all previous hits. This is a very
specialized function that helps us to merge adjacent hits into a single hit.

=over 4

=item contigs

Reference to a hash keyed on contig ID that maps each to a list of L<BasicLocation> objects for hits.

=item match

An L<Hsp> object representing a new hit.

=item RETURN

Returns TRUE if the hits were merged, else FALSE.

=back

=cut

sub MergeMatch {
    my ($self, $contigs, $match) = @_;
    # Declare the return value. This will be set to TRUE if we merge.
    my $retVal;
    my $contig = $match->sid;
    if (! exists $contigs->{$contig}) {
        # First hit on this contig. Save it as a location.
        $contigs->{$contig} = [$match->sloc];
    } else {
        # Check to see if we belong next to the previous hit.
        my $oldLoc = pop @{$contigs->{$contig}};
        my $newLoc = $match->sloc;
        if ($oldLoc->Dir eq $newLoc->Dir && ($newLoc->Left - $oldLoc->Right) < $self->{gap} ) {
            # Yes we do.
            $oldLoc->Merge($newLoc);
            push @{$contigs->{$contig}}, $oldLoc;
            $retVal = 1;
        } else {
            push @{$contigs->{$contig}}, $oldLoc, $newLoc;
        }
    }
    # Return the merge flag.
    return $retVal;
}


1;