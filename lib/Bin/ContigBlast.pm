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


package Bin::ContigBlast;

    use strict;
    use warnings;
    use BlastInterface;
    use FIG_Config;
    use File::Copy::Recursive;
    use Hsp;
    use BasicLocation;
    use FastA;


=head1 Contig Blasting Manager

This package processes blast databases for use in blasting multiple contigs from a genetic sample. The constructor takes
as input a list of FASTA files that are combined to build the blast database. If a single FASTA file is specified, it is
used as is; otherwise, a temporary file is created.

The goal is to get a list of sufficient matches for each incoming contig against the protein database.
Matches that are physically close together and target the same protein are combined. To be sufficient,
a match must have a length that is a sufficiently high proportion of the matched protein (default 60%)
and a sufficiently high e-value (default 1e-10). BLASTX is used-- the contigs are DNA, but the blast databases
are proteins.

Contigs can be processed in batches. Each batch of contigs functions as a query against the target blast database
and results are returned on a contig-by-contig basis.

The following fields are found in this object.

=over 4

=item maxE

The maximum permissible e-value for a match.

=item minLen

The minimum fraction of the protein length that must match.

=item tempDir

The name of the temporary directory to contain intermediate files.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits from the same query sequence in the same direction
that match parts of the same protein and are closer than this distance are combined.

=item blastDB

The file name of the blast database.

=item tempFlag

If TRUE, then the blast database is a temporary file and should be deleted when this object is destroyed.

=item logH

An open file handle for a log file to be used for progress messages, or C<undef> if no logging is desired.

=item stats

A statistics object for tracking statistical information.

=item batchSize

The suggested size of a batch to submit to BLAST, in base pairs.

=back

=head2 Special Methods

=head3 new

    my $blastMgr = ContigBlast->new(@fastas, \%options);

Create a new blast manager for the specified FASTA files.

=over 4

=item fastas

A list of FASTA file names. The files must be protein FASTA. Each sequence ID is the protein's feature ID, and the comment
should be the genome ID and name separated by a tab.

=item options (optional)

Reference to a hash containing zero or more of the following options.

=over 8

=item maxE

The maximum permissible e-value for a match. The default is C<1e-10>.

=item minLen

The minimum fraction of the protein length that must match.  The default is C<0.60>.

=item tempDir

The name of the temporary directory to contain intermediate files. The default is the SEEDtk temporary directory.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits from the same query sequence in the same direction
that match parts of the same protein and are closer than this distance are combined. The default is C<300>.

=item logH

Open handle to an output file for logging. The default is not to log.

=item stats

A L<Stats> object for tracking statistics. If not specified, one will be created.

=item batchSize

The suggested size of a batch to submit to BLAST, in base pairs. The default is 2000000.

=back

=back

=cut

sub new {
    my ($class, @fastas) = @_;
    # Check for an options hash.
    my $options = {};
    my $n = scalar(@fastas) - 1;
    if ($n >= 0 && ref $fastas[$n] eq 'HASH') {
        $options = pop @fastas;
    }
    # Compute the options.
    my $tempDir = $options->{tempDir} // $FIG_Config::temp;
    my $maxE = $options->{maxE} // 1e-10;
    my $minLen = $options->{minLen} // 0.60;
    my $gap = $options->{gap} // 300;
    my $logH = $options->{logH};
    my $stats = $options->{stats} // Stats->new();
    my $batchSize = $options->{batchSize} // 2000000;
    $n = scalar(@fastas);
    if ($n == 0) {
        die "There must be at least one FASTA file.";
    } else {
        for my $fasta (@fastas) {
            if (! -f $fasta) {
                die "FASTA file $fasta not found or invalid.";
            }
        }
    }
    # Now we need to create the blast DB file.
    my ($blastDB, $tempFlag);
    if (! -d $tempDir) {
        File::Copy::Recursive::pathmk($tempDir) || die "Could not create $tempDir: $!";
    }
    if ($n == 1) {
        $blastDB = $fastas[0];
        # Denote that we do not own the DB file.
        $tempFlag = 0;
        $stats->Add(fastaIn => 1);
    } else {
        $blastDB = "$tempDir/fasta$$.fa";
        open(my $oh, '>', $blastDB) || die "Could not open output fasta $blastDB: $!";
        for my $fasta (@fastas) {
            open(my $ih, '<', $fasta) || die "Could not open fasta $fasta: $!";
            while (! eof $ih) {
                my $line = <$ih>;
                print $oh $line;
            }
            $stats->Add(fastaIn => 1);
        }
        # Denote that the DB file is temporary.
        $tempFlag = 1;
    }
    # Finally, create the object.
    my $retVal = {
        maxE => $maxE,
        minLen => $minLen,
        tempDir => $tempDir,
        gap => $gap,
        blastDB => $blastDB,
        tempFlag => $tempFlag,
        logH => $logH,
        stats => $stats,
        batchSize => $batchSize,
    };
    bless $retVal, $class;
    return $retVal;
}

=head2 Public Manipulation Methods

=head3 FindMatches

    my $retHash = $blastMgr->FindMatches(\%contigs, \%retHash);

Find matches from the blast database in the specified contigs.

=over 4

=item contigs

Reference to a hash mapping each contig ID to its DNA sequence.

=item retHash (optional)

Reference to a hash used to return results, keyed on contig ID. Each contig ID will map to a 6-tuple for the best-matching protein,
consisting of (0) the protein ID, (1) a L<BasicLocation> object for the hit location, (2) the length matched,
(3) the number of match positions, (4) the genome ID, and (5) the genome name. If this parameter is not provided, a hash
will be created and returned.

=back

=cut

sub FindMatches {
    my ($self, $contigs, $retVal) = @_;
    # Insure we have a result area.
    $retVal //= {};
    # Get the gap length and the minimum length. Note that the minimum length is a fraction of the subject sequence length.
    my $gap = $self->{gap};
    my $minLen = $self->{minLen};
    my $stats = $self->{stats};
    # Create the triples for the contigs.
    my @triples = map { [$_, '', $contigs->{$_}] } keys %$contigs;
    my $inCount = scalar @triples;
    # Call the blast interface.
    $self->log("Blasting $inCount sequences.");
    my $time = time;
    my @results = BlastInterface::blastx(\@triples, $self->{blastDB},
            { outForm => 'hsp', maxE => $self->{maxE} });
    $stats->Add(blastTime => (time - $time));
    $stats->Add(blastSeqs => $inCount);
    my $hitCount = scalar @results;
    $stats->Add(blastHits => $hitCount);
    $self->log("$hitCount matches found.");
    # We now have a list of HSP matches. We need to merge them and organize them by contig. We assign a virtual
    # score to each match we will use to sort them at the end, equal to the identity count. When merging two matches
    # the counts are added. This hash maps each contig to a sub-hash of proteins to matches. A match is
    # represented as the tuple (contigLocation, proteinLength, virtualScore).
    my %contigM;
    # This tracks the total length for each protein.
    my %protL;
    # Now fill the hashes.
    for my $match (@results) {
        my $protID = $match->sid;
        my $contigID = $match->qid;
        my $contigLocation = BasicLocation->new(join('_', $contigID, $match->q1, $match->q2));
        my $protLen = $match->n_mat;
        $protL{$protID} = $match->slen;
        my ($genome, $name) = ($match->sdef =~ /^(\S+)\s+(.+)/);
        push @{$contigM{$contigID}{$protID}}, [$contigLocation, $protLen, $match->n_id, $genome, $name];
        $stats->Add(matchProcessed => 1);
    }
    # We don't need the result array any more. Everything is in hashes.
    @results = ();
    # Now we merge close matches.
    for my $contig (keys %contigM) {
        # This will hold the matches for this contig.
        my @contigHits;
        # Get the protein sub-hash.
        my $protH = $contigM{$contig};
        for my $prot (keys %$protH) {
            # Sort the matches by contig position.
            my @matches = sort { $a->[0]->Left <=> $b->[0]->Left } @{$protH->{$prot}};
            # Run through all the matches, merging them.
            my @merged;
            my $current = shift @matches;
            while (my $new = shift @matches) {
                my $merged = 0;
                if ($new->[0]->Dir eq $current->[0]->Dir) {
                    if ($new->[0]->Left - $current->[0]->Right <= $gap) {
                        # Here we can merge.
                        $stats->Add(matchMerged => 1);
                        # Merge the locations.
                        $current->[0]->Combine($new->[0]);
                        # Sum the scores.
                        $current->[2] += $new->[2];
                        # Merge the lengths.
                        $current->[1] += $new->[1];
                        $merged = 1;
                    }
                }
                if (! $merged) {
                    # Not mergable. Save the old match and keep the new one.
                    push @merged, $current;
                    $current = $new;
                    $stats->Add(matchKept => 1);
                }
            }
            # Save the last match.
            push @merged, $current;
            my $hitsLeft = scalar @merged;
            # Compute the minimum permissible match length for this protein.
            my $minProtLen = $minLen * $protL{$prot};
            # Save the acceptable match with the best score.
            my ($best) = sort { $b->[2] <=> $a->[2] } grep { $_->[1] >= $minProtLen } @merged;
            if ($best) {
                push @contigHits, [$prot, @$best];
                $stats->Add(matchDiscarded => ($hitsLeft - 1));
            } else {
                $stats->Add(matchDiscarded => $hitsLeft);
            }
        }
        # Now we have all the best hits to each hit location on this contig. If there are any, sort them by score.
        if (@contigHits) {
            $stats->Add(contigsHit => 1);
            my $best = pop @contigHits;
            while (my $new = pop @contigHits) {
                if ($new->[3] > $best->[3]) {
                    $best = $new;
                }
            }
            $retVal->{$contig} = $best;
        } else {
            $stats->Add(contigsNotHit => 1);
        }
    }
    # Return the hash of hits.
    return $retVal;
}


=head3 FindAllMatches

    my $matchHash = $blastMgr->FindAllMatches($fileName);

Find matches from the blast database in the specified contig FASTA file.

=over 4

=item fileName

The name of a DNA FASTA file for the set of contigs to BLAST.

=item RETURN

Returns a reference to a hash keyed on contig ID. Each contig ID will map to a 4-tuple for the best-matching protein,
consisting of (0) the protein ID, (1) a L<BasicLocation> object for the hit location, (2) the length matched,
and (3) the number of match positions.

=back

=cut

sub FindAllMatches {
    my ($self, $fileName) = @_;
    # This will be the return hash.
    my %retVal;
    # Get the batch size.
    my $size = $self->{batchSize};
    # Set up to read through the FASTA file.
    my $fh = FastA->new($fileName);
    # This will track the current set of contigs.
    my %contigs;
    # This will track the number of DNA letters so far.
    my $bases = 0;
    # Loop through the contigs.
    while ($fh->next) {
        my $newLen = length($fh->left);
        if ($newLen + $bases > $size || $fh->at_end) {
            # No more room. Process this batch.
            $self->FindMatches(\%contigs, \%retVal);
            # Set up for the next batch.
            %contigs = ();
            $bases = 0;
        }
        # Store this contig.
        $contigs{$fh->id} = $fh->left;
        $bases += $newLen;
    }
    # Process any residual.
    if (keys %contigs) {
        $self->FindMatches(\%contigs, \%retVal);
    }
    # Return the results.
    return \%retVal;
}


=head2 Query Methods

=head3 stats

    my $stats = $blastMgr->stats;

Return the statistics object.

=cut

sub stats {
    my ($self) = @_;
    return $self->{stats};
}


=head2 Internal Utilities

=head3 log

    $blastMgr->log($message);

Write a message to the log file.

=over 4

=item message

The message to write. A new-line will be appended.

=back

=cut

sub log {
    my ($self, $message) = @_;
    my $lh = $self->{logH};
    if ($lh) {
        print $lh "$message\n";
    }
}

1;