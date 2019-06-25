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


package KmerFramerFiles;

    use strict;
    use base qw(KmerFramer);

=head1 Write Coding Frame Kmers to a File

This is a subclass of L<KmerFramer>.  Instead of storing the kmers in a hash in memory, they are written
to the standard output in tab-delimited format.  Each record of the standard output will contain two columns
(0) the kmer and (1) the frame ID.

=head2 Internal Utilities

=head3 _StoreKmer

    $kmerFramer->_StoreKmer($kmer, $frame);

Output the kmer record for the specified frame.

=over 4

=item kmer

Kmer sequence in the forward direction.

=item frame

Frame identifier for the kmer.

=back

=cut

sub _StoreKmer {
    my ($self, $kmer, $frame) = @_;
    my $stats = $self->{stats};
    # Get the reverse complement.
    my $kmerR = SeedUtils::rev_comp($kmer);
    for my $kmerX ($kmer, $kmerR) {
        print join("\t", $kmerX, $frame) . "\n";
        $frame =~ tr/+-/-+/;
    }
}


1;

