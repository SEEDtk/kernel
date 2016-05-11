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


package Bin::Kmer;

    use strict;
    use warnings;
    use base qw(KmerDb);

=head1 Assign Bins via Kmer Discrimination

This module is used to create a kmer discrimination database for binning. All of the bin's reference
genomes are fed in, and kmers created from the attendant protein sequences. The L</accumulate_hits>
method is overridden to just count kmers that occur in only one bin. A separate method uses the kmers to 
determine which bin should contain any incoming contig. Note that the incoming contigs contain DNA while 
the kmers are proteins, so translation is necessary.

In the parlance of the base class L<KmerDb>, the groups are bins, and the sources are incoming contigs.

This object has the following fields, in addition to those in the base object L<KmerDb>.

=over 4

=item binStrength

Number of matches required to place a sequence into a bin. The default is C<10>.

=back

=head2 Special Methods

=head3 new

    my $kmerdb = Bin::Kmer->new(%options);

Create a new, blank KMER database. The following options are supported.

=over 4

=item kmerSize

The size of a kmer. The default is C<10>.

=item binStrength

Number of matches required to place a sequence into a bin. The default is C<10>.

=item json

The name of a JSON file (or an open file handle) containing the group and kmer hashes. This file is created by
the L</Save> method.

=back

=cut

sub new {
    my ($class, %options) = @_;
    # Create the object.
    my $retVal = $class->SUPER::new(%options);
    # Store the bin strength.
    $retVal->{binStrength} = $options{binStrength} // 10;
    # Return the object.
    return $retVal;
}

=head2 Public Methods

=head3 AddGenome

    $kmerdb->AddGenome($gto, $bin);

Add a reference genome for a bin. Kmers will be generated for all of the genome's proteins.

=over 4

=item gto

L<GenomeTypeObject> for a bin's reference genome. (Some bins may have more than one.)

=item bin

ID of the bin.

=back

=cut

sub AddGenome {
    my ($self, $gto, $bin) = @_;
    # Loop through the proteins for the GTO, converting them to kmers.
    my $flist = $gto->{features};
    for my $fdata (@$flist) {
        # Insure we have a protein for this feature.
        my $sequence = $fdata->{protein_translation};
        if ($sequence) {
            # We do, get kmers from it.
            $self->AddSequence($bin, $sequence);
        }
    }
}

=head3 Finalize

    $kmerdb->Finalize();

Clean up the kmer hash to remove kmers found in more than one group. This prepares the database
for queries, and overrides the base class method. Because computing discriminators is very
expensive, we simply mark the bin as finalized and override accumulate_hits to only count kmers in
a single bin.

=cut

sub Finalize {
    my ($self) = @_;
    $self->{finalized} = 1;
}



=head3 ComputeBin

    my $groupID = $self->ComputeBin($sequence, %options);

Compute the bin ID for a specified sequence.

=over 4

=item sequence

Sequence to be examined.

=item options

Zero or more of the following options.

=over 8

=item translate

If specified, the genetic code for translating the sequence. The incoming sequence is presumed to be DNA and the
kmers are presumed to be proteins. The incoming sequence will be translated before it is scanned for kmer matches.

=back

=item RETURN

Returns the ID of the bin to which the sequence belongs, or C<undef> if the sequence does not match a particular
bin.

=back

=cut

sub ComputeBin {
    my ($self, $sequence, %options) = @_;
    # This hash will count the number of discriminating kmers found for each bin in the incoming sequence.
    my %binHash;
    # Compute the number of hits.
    $self->count_hits($sequence, \%binHash, $options{translate});
    # Find the bin with the most hits greater than or equal to the bin strength.
    my ($bestBin, $bestCount) = ('', $self->{binStrength} - 1);
    my $errorBin = '';
    for my $bin (keys %binHash) {
        if ($binHash{$bin} > $bestCount) {
            ($bestBin, $bestCount) = ($bin, $binHash{$bin})
        } elsif ($binHash{$bin} == $bestCount) {
            $errorBin = $bin;
            $bestBin = $bin;
        }
    }
    # This will be the return value.
    my $retVal;
    # Check for a tie at the top.
    if ($bestBin ne $errorBin) {
        # Here we did NOT have a tie for the best bin. We take that as confirmation we've made the right choice.
        $retVal = $bestBin;
    }
    # Return the bin found.
    return $retVal;
}

=head2 Internal Methods

=head3 accumulate_hits

    $kmerdb->accumulate_hits($sequence, \%counts);

Accumulate the kmer hits in a sequence. The sequence must already be translated if translation is
needed. This is an override to the base-class method that only counts a hit if it has just a single
target bin.

=over 4

=item sequence

Sequence to examine for kmers.

=item counts

Reference to a hash mapping group IDs to hit counts.

=back

=cut

sub accumulate_hits {
    my ($self, $sequence, $counts) = @_;
    # Get the kmer length and the hash.
    my $kmerSize = $self->{kmerSize};
    my $kmerHash = $self->{kmerHash};
    # Loop through the kmers.
    my $n = length($sequence) - $kmerSize;
    for (my $i = 0; $i <= $n; $i++) {
        my $kmer = substr($sequence, $i, $kmerSize);
        my $groups = $kmerHash->{$kmer};
        if ($groups) {
            my ($group, @others) = keys %$groups;
            if (! @others) {
                $counts->{$group}++;
            }
        }
    }
}


1;