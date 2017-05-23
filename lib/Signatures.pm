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


package Signatures;

    use strict;
    use warnings;
    use SeedUtils;

=head1 Manage a Protein Signatures Database

This package manages an in-memory protein signatures database. The database is a hash of protein signature
kmers matching each to a given group ID (usually a representative genome ID). The binning algorithm takes as
input a sequence reader object-- either L<FastA> or L<FastQ>-- and counts the signatures. The signature
group with the highest score gets custody of the sequences.

This package is a virtual class. The constructor must specify a subclass that has a scoring algorithm defined
(implementation of the L</score> function).

This object contains the following fields.

=over 4

=item sigsH

Reference to a hash mapping each signature kmer to its group ID.

=item groupH

Reference to a hash mapping each group ID to the number of signatures for it.

=item kmerSize

The kmer size for the signatures.

=item xlate

A translation table (from L<SeedUtils/genetic_code>)to use for protein translation.

=item sigCount

Total number of signatures.

=item groupCount

Total number of groups.

=item minScore

Minimum permissible score for a sequence to be binned.

=item positions

The total number of kmer positions processed for the current sequence.

=item stats

A L<Stats> object for keeping statistics.

=back

=head2 Special Methods

=head3 new

    my $sigDB = Signatures->new($inFile, $min, $stats, %options);

Create a new signatures database.

=over 4

=item inFile

Name of the file containing the signatures, in tab-delimited format, each record consisting of (0) a protein signature kmer, and
(1) a group ID. An open file handle can also be specified.

=item min

Minimum permissible score for a sequence to be binned.

=item stats

A L<Stats> object for keeping statistics.

=item options

A hash of options, including zero or more of the following keys.

=over 8

=item geneticCode

The protein translation code. The default is C<11>.

=back

=back

=cut

sub new {
    my ($class, $inFile, $min, $stats, %options) = @_;
    # Get the options.
    my $geneticCode = $options{geneticCode} // 11;
    # Initialize the hashes.
    my (%sigsH, %groupH);
    # Initialize the counters.
    my ($groupCount, $sigCount) = (0, 0);
    # Initialize the kmer size holding.
    my $kmerSize;
    # Open the file for input.
    my $ih;
    if (ref $inFile eq 'GLOB') {
        $ih = $inFile;
    } else {
        open($ih, "<$inFile") || die "Could not open signatures file $inFile: $!";
    }
    # Loop through the file, filling the hashes.
    while (! eof $ih) {
        my $line = <$ih>;
        $line =~ s/[\r\n]+$//;
        my ($kmer, $group) = split /\t/, $line;
        $sigsH{$kmer} = $group;
        $sigCount++;
        $stats->Add(sigKmersRead => 1);
        if (exists $groupH{$group}) {
            $groupH{$group} = 1;
            $groupCount++;
            $stats->Add(sigGroupsRead => 1);
        } else {
            $groupH{$group}++;
        }
        my $len = length $kmer;
        if (! defined $kmerSize) {
            $kmerSize = $len;
        } elsif ($kmerSize != $len) {
            die "Invalid kmer \"$kmer\" found (length $len instead of $kmerSize.";
        }
    }
    if (! defined $kmerSize) {
        die "Kmer signatures file was empty or invalid.";
    }
    # Compute the genetic code translation table.
    my $xlate = SeedUtils::genetic_code($geneticCode);
    # Create the object.
    my $retVal = {
        sigsH => \%sigsH,
        groupH => \%groupH,
        sigCount => $sigCount,
        groupCount => $groupCount,
        xlate => $xlate,
        kmerSize => $kmerSize,
        minScore => $min,
        stats => $stats,
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Query Methods

=head3 ComputeBin

    my @binIDs = $sigDB->ComputeBin($seqObject);

Compute the bin to contain the specified sequence object. We count the kmers for each group, then
ask the scoring algorithm to compute the scores. The group with the highest score greater than the
minimum is the group chosen. If there is a tie, all of the tied bins are returned.

=over 4

=item seqObject

A L<FastA> or L<FastQ> object containing the sequences.

=item RETURN

Returns a list of viable bin IDs.

=back

=cut

sub ComputeBin {
    my ($self, $seqObject) = @_;
    # Get the list of sequences to process.
    my @seqs = $seqObject->seqs;
    # Get the statistics object.
    my $stats = $self->{stats};
    # Initialize the length counter for this computation.
    $self->{positions} = 0;
    # This will count our hits.
    my %hits;
    # Loop through the sequences, counting hits. Note we do both the sequence and its
    # reverse complement.
    for my $seq (@seqs) {
        $stats->Add(dnaSequences => 1);
        $self->_CountDnaHits(\$seq, \%hits);
        SeedUtils::rev_comp(\$seq);
        $self->_CountDnaHits(\$seq, \%hits);
    }
    # Get the group counts hash.
    my $groupH = $self->{groupH};
    # Get the position count.
    my $positions = $self->{positions};
    # Use the counts to compute scores.
    my @retVal;
    my $best = $self->{minScore};
    for my $group (keys %hits) {
        my $score = $self->score($hits{$group}, $groupH->{$group}, $positions);
        if ($score > $best) {
            @retVal = ($group);
        } elsif ($score == $best) {
            push @retVal, $group;
        }
    }
    my $found = scalar(@retVal);
    $stats->Add("binsFound$found" => 1);
    # Return the bins found.
    return @retVal;
}

=head2 Internal Methods

=head3 _CountDnaHits

    $sigsObject->_CountDnaHits(\$seq, \%hits);

Count the signature hits in a DNA sequence. There is no need to reverse-complement, but we do need to translate
to proteins.

=over 4

=item seq

A reference to the DNA sequence to process.

=item hits

A reference to a hash where the hit counts are being stored.

=back

=cut

sub _CountDnaHits {
    my ($self, $seq, $hits) = @_;
    # Get the translation matrix.
    my $xlate = $self->{xlate};
    # Loop through the frames.
    for my $f (0, 1, 2) {
        my $protSeq = SeedUtils::translate($seq, $f, $xlate);
        $self->_CountProtHits(\$protSeq, $hits);
    }
}

=head3 _CountProtHits

    $self->_CountProtHits(\$seq, \%hits);

Count the hits in the specified protein sequence.

=over 4

=item seq

Reference to a protein sequence.

=item hits

Reference to the hash in which the hit counts should be placed.

=back

=cut

sub _CountProtHits {
    my ($self, $seq, $hits) = @_;
    # Get the kmer size.
    my $kmerSize = $self->{kmerSize};
    # Get the statistics object.
    my $stats = $self->{stats};
    # Get the kmer hash.
    my $sigsH = $self->{sigsH};
    # Loop through the sequence.
    my $len = length($$seq);
    my $last = $len - $kmerSize;
    for (my $i = 0; $i < $last; $i++) {
        # Get this kmer and count it.
        my $kmer = substr($$seq, $i, $kmerSize);
        $stats->Add(kmersChecked => 1);
        $self->{positions}++;
        # Is it a signature?
        my $group = $sigsH->{$kmer};
        if ($group) {
            # Yes. Count the hit.
            $stats->Add(kmersCounted => 1);
            $hits->{$group}++;
        }
    }
}

=head2 Virtual Methods

=head3 score

    my $score = $sigsObject->score($hits, $groupCount, $positions);

Compute the score for this group given the specified number of hits.

=over 4

=back

=cut

sub score {
    my ($self, $hits, $groupCount, $positions) = @_;
    my $retVal;
    die "Pure virtual method \"score\" called.";
    return $retVal;
}


1;