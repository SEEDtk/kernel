#
# Copyright (c) 2003-2019 University of Chicago and Fellowship
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


package KmerFramer;

    use strict;
    use warnings;
    use Stats;
    use P3DataAPI;
    use P3Utils;
    use SeedUtils;
    use Statistics::Descriptive;
    use BasicLocation;
    use Shrub::Contigs;
    use POSIX qw(ceil floor);

=head1 Analyze Sequences to Classify Kmers by Frame

This object creates a kmer database that counts the number of times each kmer occurs within a specific coding frame.
The coding frames are C<-3>, C<-2>, C<-1>, C<0>, C<+1>, C<+2>, and C<+3>.  C<0> indicates a non-coding region.  The
others indicate the strand (is the DNA on the same or opposite strand as the protein coding) and whether the DNA begins
in the same codon frame (C<3>), or is off by 1 or 2.

The fields in this object are as follows.

=over 4

=item p3

A L<P3DataAPI> object for accessing the PATRIC database.

=item debug

If specified, an open file handle on which to write progress messages.

=item stats

A L<Stats> object used to track activity.

=item K

The DNA kmer size.

=item kHash

Reference to a hash mapping each kmer to a 7-tuple containing the counts, in the order C<-3>, C<-2>, C<-1>, C<0>, C<+1>,
C<+2>, and C<+3>.

=item gHash

Reference to a hash keyed on the IDs of all the genomes processed for kmers.

=back

=head2 Special Methods

=head3 new

    my $kmerFramer = KmerFramer->new(%options);

Create a new, empty kmer region analysis database.

=over 4

=item options

A hash containing zero or more of the following keys.

=over 8

=item p3

A L<P3DataAPI> object for accessing the PATRIC database.  If none is provided, one will be created.
Alternatively, a L<Shrub> object for accessing the Shrub database.

=item K

The DNA kmer size.  The default is C<15>.

=item debug

An open file handle on which to write progress messages.  If none is specified, no progress messages will be written.

=item stats

A L<Stats> object for tracking activity.  If none is provided, one will be created.

=item saved

If specified, the name of a file containing a saved copy of the database (created using L</Save>).  If this option is
specified, the C<K> value is ignored.

=back

=back

=cut

sub new {
    my ($class, %options) = @_;
    # Get the options.
    my $debug = $options{debug};
    my $p3 = $options{p3} // P3DataAPI->new();
    my $stats = $options{stats} // Stats->new();
    my $K = $options{K} // 15;
    # This will be the return value.
    my $retVal;
    # Are we loading?
    if ($options{saved}) {
        # Yes.  Create the object from the file.
        $retVal = _Load($options{saved});
    } else {
        # No.  Create a blank object.
        my (%kHash, %gHash);
        $retVal = { kHash => \%kHash, gHash => \%gHash, K => $K };
    }
    # Update the run-time options.
    $retVal->{debug} = $debug;
    $retVal->{p3} = $p3;
    $retVal->{stats} = $stats;
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 INDICES

This is a constant hash that maps each frame label to its position in the kmer hash list.

=cut

use constant INDICES => { '-3' => 0, '-2' => 1, '-1' => 2, '0' => 3, '+1' => 4, '+2' => 5, '+3' => 6 };

=head3 LABELS

This is a constant list that maps each position in the kmer hash list to its frame label.

=cut

use constant LABELS => ['-3', '-2', '-1', '0', '+1', '+2', '+3'];

=head2 Access Methods

=head3 stats

    my $stats = $kmerFramer->stats;

Return the internal L<Stats> object.

=cut

sub stats {
    my ($self) = @_;
    return $self->{stats};
}

=head3 debug

    my $debug = $kmerFramer->debug;

Return the internal debug log handle.

=cut

sub debug {
    my ($self) = @_;
    return $self->{debug};
}


=head2 Public Manipulation Methods

=head3 CountKmers

    $kmerFramer->CountKmers($sequence, $seqMap);

This method counts the kmers in a sequence using a coding region map output by L</SequenceMap>.

=over 4

=item sequence

The sequence whose kmers are to be counted.

=item seqMap

A reference to the list of coding regions in the sequence.  Each element of the list is a 3-tuple containing (0) the
start location (1-based), (1) the end location (1-based), and (2) the strand (C<+>, C<->, or C<0> for invalid).  The
regions are non-overlapping, and sorted by start location.

=back

=cut

sub CountKmers {
    my ($self, $sequence, $seqMap) = @_;
    my $K = $self->{K};
    my $stats = $self->stats;
    # Get the sequence length.
    my $seqLen = length $sequence;
    # This map entry will be used when we run off the end.
    my $trailer = [$seqLen+1, $seqLen+1, '0'];
    # Pick off the first map entry.
    my $mEntry = $seqMap->[0] // $trailer;
    my $mPos = 1;
    # Compute the location of the last kmer.  Note we are using 1-based locations.
    my $last = $seqLen - $K + 1;
    # Loop through the sequence.
    for (my $i = 1; $i <= $last; $i++) {
        # Compute this kmer's end location.
        my $end = $i + $K - 1;
        # Get the sequence.
        my $kmer = uc substr($sequence, $i-1, $K);
        # Loop until we find what this kmer is.
        my $found;
        while (! $found) {
            # Is the whole kmer in front of the current region?
            if ($end < $mEntry->[0]) {
                # Store it as non-coding.
                $self->_StoreKmer($kmer, '0');
                $found = 1;
            } elsif ($i < $mEntry->[0]) {
                # Here we overlap, so we ignore this one.
                $stats->Add(boundaryKmer => 1);
                $found = 1;
            } elsif ($end <= $mEntry->[1]) {
                # Here we are inside the coding region.  Is the coding region valid?
                my $strand = $mEntry->[2];
                if ($strand eq '0') {
                    # No.  Ignore this kmer.
                    $stats->Add(ignoreKmer => 1);
                } else {
                    # Compute the frame.
                    my $frm = $strand . ((($i - $mEntry) % 3) || 3);
                    $self->_StoreKmer($kmer, $frm);
                }
                $found = 1;
            } elsif ($i <= $mEntry->[1]) {
                # Here we overlap again, so we ignore.
                $stats->Add(boundaryKmer => 1);
                $found = 1;
            } else {
                # Here we are past the end of this region.  Pop the next entry off the sequence map.
                $mEntry = $seqMap->[$mPos] // $trailer;
                $mPos++;
            }
        }
    }
}

=head3 Save

    $kmerFramer->Save($outFile);

Store this kmer database to a file in JSON format.

=over 4

=item outFile

The name of the output file, or an open file handle to which the database should be written.

=back

=cut

sub Save {
    my ($self, $outFile) = @_;
    my $oh;
    if (ref $outFile eq 'GLOB') {
        $oh = $outFile;
    } else {
        open($oh, '>', $outFile) || die "Could not open JSON output file for kmers: $!";
        print $oh $self->{K} . "\n";
        my $kHash = $self->{kHash};
        for my $kmer (keys %$kHash) {
            print $oh join("\t", $kmer, @{$kHash->{$kmer}}) ."\n";
        }
        my $gHash = $self->{gHash};
        for my $genome (keys %$gHash) {
            print $oh "$genome\n";
        }
    }
}

=head3 Record

    $kmerFramer->Record($genome);

Denote that a genome has been processed.

=over 4

=item genome

ID of the genome that has been processed.

=back

=cut

sub Record {
    my ($self, $genome) = @_;
    $self->{gHash}{$genome} = 1;
}



=head2 Database Methods

=head3 SequenceMap

    my $seqMap = $kmerFramer->SequenceMap($genomeID);

Create a sequence map for the specified genome.  The sequence map is reference to a hash of all the non-overlapping
coding regions in the genomes, keyed by sequence ID and sorted by starting offset within sequence ID.  The information is
read from PATRIC, and can be used to compute the frame of a kmer in the genome.  In practice, we will run through the
genome's sequences in alphabetical order by sequence ID, and then from left to right on the sequence, and we will
run through the sequence maps in parallel.  Note that regions containing overlapping features, or frame-shifted features,
are not counted.  These are marked appropriately.

=over 4

=item genomeID

The ID of the genome whose sequences are to be analyzed.

=item RETURN

Returns a reference to a hash mapping each sequence ID to a list of 3-tuples.  Each 3-tuple consists of (0) a starting
offset (0-based), an ending offset, and a strand indicator-- C<+> if there is a protein on the + strand, C<-> if there
is a protein on the - strand, and C<0> if the region is invalid and should not be counted.  Note that all the offsets
are 1-based.

=back

=cut

sub SequenceMap {
    my ($self, $genomeID) = @_;
    my $p3 = $self->{p3};
    my $stats = $self->stats;
    my $debug = $self->debug;
    # Get the genome's CDS features.
    print $debug "Reading coding features for $genomeID.\n" if $debug;
    my $cdsList;
    if (ref $p3 eq 'P3DataAPI') {
        # For PATRIC we get all the features en masse.
        $cdsList = P3Utils::get_data($p3, feature => [['eq', 'genome_id', $genomeID], ['eq', 'feature_type', 'CDS']],
                ['patric_id', 'sequence_id', 'strand', 'start', 'end', 'segments']);
    } else {
        $cdsList = [];
        # For Shrub, we loop through the features and must compute the global strand, start, and end along with the
        # segments.
        my @fids = $p3->GetFlat('Genome2Feature', 'Genome2Feature(from-link) = ? AND Genome2Feature(to-link) LIKE ?',
                [$genomeID, "fig|$genomeID.peg.%"], 'to-link');
        for my $fid (@fids) {
            my @locs = $p3->fid_locs($fid);
            my ($start, $end) = BasicLocation::ListBounds(@locs);
            my @segments = map { [$_->Left, $_->Right] } @locs;
            my $sequenceID = $locs[0]->Contig;
            my $strand = $locs[0]->Dir;
            push @$cdsList, [$fid, $sequenceID, $strand, $start, $end, \@segments];
        }
    }
    my $fCount = scalar @$cdsList;
    print $debug "Sorting $fCount results.\n" if $debug;
    # Separate them by sequence.
    my %seqHash;
    for my $cdsItem (@$cdsList) {
        my ($id, $seqID, $strand, $start, $end, $segments) = @$cdsItem;
        # If this is a multi-segment feature, it counts as the zero strand.
        if ($segments && @$segments > 1) {
            $strand = '0';
            $stats->Add(frameShiftFound => 1);
        }
        push @{$seqHash{$seqID}}, [$start, $end, $strand];
    }
    # Release the master list.
    undef $cdsList;
    # This will be the return hash.
    my %retVal;
    # Sort each sequence.
    for my $seqID (keys %seqHash) {
        print $debug "Processing features in $seqID.\n" if $debug;
        my $seqList = $seqHash{$seqID};
        $seqList = [sort { $b->[0] <=> $a->[0] } @$seqList];
        # Now we have the sequence's entries sorted.  We need to transfer them to the output, separating out the
        # overlaps.  We will be popping and pushing items onto the active stack, and putting finished items
        # into the result queue.
        my $firstItem = pop @$seqList;
        my ($start, $end, $strand) = @$firstItem;
        # Loop through the items, transferring them to the output list.
        my @retList;
        while (my $seqItem = pop @$seqList) {
            my ($start1, $end1, $strand1) = @$seqItem;
            if ($start1 <= $end) {
                $stats->Add(overlap => 1);
                # Here we have an overlap.
                if ($start1 > $start) {
                    # Hack off the non-overlapping part.
                    push @retList, [$start, $start1-1, $strand];
                    $stats->Add(overlapPrefixOut => 1);
                }
                # Now we know everything about the region before $start1.  $start1 is the beginning of an overlap region.
                $start = $start1;
                my $strand0 = $strand;
                $strand = '0';
                # The overlap region ends at end1 or end, whichever is smaller.  Everything after that is coding on the
                # appropriate item's strand.
                if ($end1 < $end) {
                    push @$seqList, [$end1+1, $end, $strand1];
                    $end = $end1;
                } elsif ($end1 > $end) {
                    push @$seqList, [$end+1, $end1, $strand0];
                }
            } else {
                # Here the regions are non-overlapping.  Output the old one and save the new one.
                push @retList, [$start, $end, $strand];
                ($start, $end, $strand) = ($start1, $end1, $strand1);
                $stats->Add(regionOut => 1);
            }
        }
        push @retList, [$start, $end, $strand];
        $retVal{$seqID} = \@retList;
    }
    return \%retVal;
}

=head3 GenomeNames

    my $nameResults = $kmerFramer->GenomeNames(\@genomeIDs);

Retrieve a batch of genome names from whichever database is attached.

=over 4

=item genomeIDs

Reference to a list of genome IDs.

=item RETURN

Returns a reference to a list of 2-tuples, each consisting of (0) a genome ID and (1) the genome name.

=back

=cut

sub GenomeNames {
    my ($self, $genomeIDs) = @_;
    my $p3 = $self->{p3};
    # This will be the return value.
    my $retVal;
    if (ref $p3 eq 'P3DataAPI') {
        # Here we have PATRIC, and we ask for all the genomes at once.
        my @couples = map { [$_, [$_]] } @$genomeIDs;
        $retVal = P3Utils::get_data_batch($p3, genome => [], ['genome_name'], \@couples);
    } else {
        # Here we have a Shrub database, and we rip through the genomes one at a time.
        $retVal = [];
        for my $genome (@$genomeIDs) {
            my ($name) = $p3->GetFlat('Genome', 'Genome(id) = ?', [$genome], 'name');
            if ($name) {
                push @$retVal, [$genome, $name];
            }
        }
    }
    return $retVal;
}

=head3 SequenceList

=cut

sub SequenceList {
    my ($self, $genomeID) = @_;
    my $p3 = $self->{p3};
    my $retVal;
    if (ref $p3 eq 'P3DataAPI') {
        # Here we can ask for the sequences from PATRIC.
        $retVal = P3Utils::get_data($p3, contig => [['eq', 'genome_id', $genomeID]], ['sequence_id', 'sequence']);
    } else {
        # For Shrub, we use a Contigs object.
        my $contigs = Shrub::Contigs->new($p3, $genomeID);
        $retVal = [map { [$_->[0], $_->[2]] } $contigs->tuples];
    }
    return $retVal;
}

=head2 Query Methods

=head3 frac

    my ($hFrac, $frame) = $kmerFramer->_frac($kmer);

For a specified kmer, return the most likely frame identifier and the fraction of occurrences
in that frame.  If the kmer does not exist, the frame will be undefined and the fraction will
be zero.

=over 4

=item kmer

The kmer to examine.

=item RETURN

Returns a two-element list consisting of (0) the probability of the kmer being in a particular
frame and (1) the ID of the frame.

=back

=cut

sub frac {
    my ($self, $kmer) = @_;
    # Get the kmer's count vector.
    my $vector = $self->{kHash}{$kmer} // [0, 0, 0, 0, 0, 0, 0];
    # These will be the return values.
    my ($hFrac, $frame) = (0);
    # This is the total for the vector.
    my $tot = 0;
    for (my $i = 0; $i < 7; $i++) {
        my $v = $vector->[$i];
        if ($v > $hFrac) {
            $frame = LABELS->[$i];
            $hFrac = $v;
        }
        $tot += $v;
    }
    # Convert the best result to a fraction.
    if ($tot > 0) {
        $hFrac /= $tot;
    }
    # Return the results.
    return ($hFrac, $frame);
}

=head3 Metrics

    my ($mean, $sdev, $count) = $kmerFramer->Metrics();

Compute the mean, standard deviation, and count of the probability values of each kmer.  The
probability is the fraction of times that the kmer is in its most likely frame.  This method
returns the mean and standard deviation over the whole kmer database.  The count will
be the total number of unique kmers.

=cut

sub Metrics {
    my ($self) = @_;
    my $kHash = $self->{kHash};
    # Create a statistical calculator.
    my $calc = Statistics::Descriptive::Sparse->new();
    for my $kmer (keys %$kHash) {
        # Get the probability
        my ($frac) = $self->frac($kmer);
        $calc->add_data($frac);
    }
    return ($calc->mean(), $calc->standard_deviation(), $calc->count());
}

=head3 Distribution

    my $distH = $kmerFramer->Distribution($mean, $sdev);

Given a mean and standard deviation, compute how many kmers are in each z-score block.
The z-score is in this case computed as the probability value less the mean divided
by the standard deviation and rounded away from 0.  Thus, values equal to the exact mean return
as 0, values less than one deviation away return as 1 or -1, and so forth.  The result gives
us a crude estimate of the shape of the distribution curve.

=over 4

=item mean

The mean highest-fraction value.

=item sdev

The standard deviation of the highest-fraction values.

=item RETURN

Returns a reference to a hash mapping each z-score bracket to its kmer count.

=back

=cut

sub Distribution {
    my ($self, $mean, $sdev) = @_;
    my $debug = $self->{debug};
    # This will be the return hash.
    my %retVal;
    # Loop through the kmers.
    my $kHash = $self->{kHash};
    for my $kmer (keys %$kHash) {
        my ($frac) = $self->frac($kmer);
        my $z = ($frac - $mean) / $sdev;
        my $bracket = ($z < 0 ? floor($z) : ceil($z));
        $retVal{$bracket}++;
    }
    # Return the hash.
    return \%retVal;
}

=head3 gCheck

    my $processed = $kmerFramer->gCheck($genome);

Return TRUE if the specified genome has already been processed, else FALSE.

=over 4

=item genome

ID of the genome to check.

=item RETURN

Returns TRUE if the genome is already processed, else FALSE.

=back

=cut

sub gCheck {
    my ($self, $genome) = @_;
    return ($self->{gHash}{$genome} ? 1 : 0);
}


=head2 Internal Utilities

=head3 _Load

    my $kmerFramer = KmerFramer::_Load($file);

Load this object from a file saved using L</Save>.

=over 4

=item file

The name of the file containing the saved kmer database.

=item RETURN

Returns an unblessed object ready to use.

=back

=cut

sub _Load {
    my ($file) = @_;
    # First read the kmer size.
    open(my $ih, '<', $file) || die "Could not open input $file: $!";
    my $K = <$ih>;
    chomp $K;
    # Now read the hashes.
    my (%kHash, %gHash);
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\d+\.\d+)/) {
            $gHash{$1} = 1;
        } else {
            chomp $line;
            my ($kmer, @counts) = split /\t/, $line;
            $kHash{$kmer} = \@counts;
        }
    }
    # Return the object.
    my $retVal = { K => $K, kHash => \%kHash, gHash => \%gHash };
    return $retVal;
}

=head3 _StoreKmer

    $kmerFramer->_StoreKmer($kmer, $frame);

Update the kmer hash for the specified kmer and the specified frame.

=over 4

=item kmer

Kmer sequence in the forward direction.

=item frame

Frame identifier for the kmer.

=back

=cut

sub _StoreKmer {
    my ($self, $kmer, $frame) = @_;
    my $kHash = $self->{kHash};
    my $stats = $self->{stats};
    # Compute the frame index.
    my $fidx = INDICES->{$frame};
    # Get the reverse complement.
    my $kmerR = SeedUtils::rev_comp($kmer);
    for my $kmerX ($kmer, $kmerR) {
        # If this is a new kmer, create an empty array for it.
        if (! $kHash->{$kmerX}) {
            $kHash->{$kmerX} = [0, 0, 0, 0, 0, 0, 0];
        }
        # Update the indicated frame.
        $kHash->{$kmerX}[$fidx]++;
        # Invert the frame for the next version.
        $fidx = 6 - $fidx;
    }
}

1;


