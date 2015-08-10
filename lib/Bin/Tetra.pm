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


package Bin::Tetra;

    use strict;
    use warnings;

=head1 Create a TetraNucleotide Vector for a Contig

This object takes DNA sequences as input and outputs a tetranucleotide vector. The DNA sequence
can be in the form of a string or input from a FASTA file.

The object contains the following fields.

=over 4

=item tetraMap

Reference to a hash that maps a tetranucleotide string (upper case) to a list of output vector positions.

=item buckets

The number of positions used in the output vector.

=back

There are three schemes for mapping tetranucleotides.

=over 4

=item raw

There are 256 possibilities. Each four-base combination maps to itself and its reverse compliment.

=item fancy

There are 256 possibilities. Each four-base combination maps to itself and its reverse compliment unless the reverse
compliment is the same, in which case it only maps to itself.

=item dual

There are 136 possibilities. Each four-base combination and its reverse compliment map to the same position.

=back

=head2 Special Methods

=head3 new

    my $tetra = Bin::Tetra->new($scheme);

Create a new tetranucleotide processing object.

=over 4

=item scheme

The scheme to use for the tetranucleotide mapping-- C<raw>, C<fancy>, or C<dual>. The default is C<dual>.

=back

=cut

sub new {
    my ($class, $scheme) = @_;
    # Default the scheme.
    $scheme //= 'dual';
    # This tracks the available vector positions.
    my $next = 0;
    # This will be our mapping hash.
    my %tetraMap;
    # Note that the nucleotide at position X is the complement of the one at 3 - X.
    my @letters = ('A','G', 'C', 'T');
    for (my $num = 0; $num < 256; $num++) {
        my $norm = "";
        my $comp = "";
        my $numWork = $num;
        for (my $i = 0; $i < 4; $i++) {
            $norm .= $letters[$numWork & 3];
            $comp = $letters[3 - ($numWork & 3)] . $comp;
            $numWork >>= 2;
        }
        # Do we already have this pair?
        if (! exists $tetraMap{$norm}) {
            # No. Process it according to the scheme.
            my @buckets;
            push @buckets, $next++;
            if ($scheme ne 'dual') {
                if ($comp ne $norm) {
                    push @buckets, $next++;
                } elsif ($scheme eq 'raw') {
                    push @buckets, $next;
                }
            }
            $tetraMap{$norm} = \@buckets;
            $tetraMap{$comp} = \@buckets;
        }
    }
    # Create  the object.
    my $retVal = {
        tetraMap => \%tetraMap,
        buckets => $next
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Public Manipulation Methods

=head3 ProcessString

    my $vector = $tetra->ProcessString($dna);

Compute the tetranucleotide vector for a DNA string.

=over 4

=item dna

A string of DNA letters.

=item RETURN

Returns a reference to a vector of the tetranucleotide counts for the DNA string.

=back

=cut

sub ProcessString {
    my ($self, $dna) = @_;
    # Get the bucket map.
    my $tetraMap = $self->{tetraMap};
    # Create the return vector,
    my @retVal = (0) x $self->{buckets};
    # Loop through the DNA.
    my $n = length($dna) - 3;
    for (my $i = 0; $i < $n; $i++) {
        my $sub4 = substr($dna, $i, 4);
        my $buckets = $tetraMap->{uc $sub4} // [];
        for my $bucket (@$buckets) {
            $retVal[$bucket]++;
        }
    }
    # Return the vector.
    return \@retVal;
}


=head3 ProcessFasta

    my $contigHash = $tetra->ProcessFasta($ih);

Read all the contigs in a FASTA input file and return a hash mapping each contig ID to its
tetranucleotide vector.

=over 4

=item ih

Open input handle for the FASTA file, or the name of the FASTA file.

=item RETURN

Returns a reference to a hash mapping each input contig ID to a vector of tetranucleotide counts.

=back

=cut

sub ProcessFasta {
    my ($self, $ih) = @_;
    # Insure we have an open file handle.
    if (ref $ih ne 'GLOB') {
        open(my $fh, "<", $ih) || die "Could not open FASTA input file $ih: $!";
        $ih = $fh;
    }
    # This will be our return hash.
    my %retVal;
    # Loop through the input file.
    my @dna;
    my $contigID;
    while (! eof $ih) {
        my $line = <$ih>;
        # Is this a header line?
        if ($line =~ /^>(\S+)/) {
            # Yes. Save the new contig ID.
            my $newContig = $1;
            # Process the current DNA, if any.
            if (@dna) {
                $retVal{$contigID} = $self->ProcessString(join("", @dna));
            }
            # Set up for the next contig.
            $contigID = $newContig;
            @dna = ();
        } else {
            # Not a header line. Save the DNA.
            chomp $line;
            push @dna, $line;
        }
    }
    # If there is leftover DNA, process it.
    if (@dna) {
        $retVal{$contigID} = $self->ProcessString(join("", @dna));
    }
    # Return the hash of contigs to vectors.
    return \%retVal;
}

1;