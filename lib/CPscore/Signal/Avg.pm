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


package CPscore::Signal::Avg;

    use strict;
    use warnings;
    use base qw(CPscore::Signal);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on
similarity signal vectors. In this case, the similarity vector contains the mean of the percent identity times the length
of all the BLAST matches to each reference genome. The resulting number is an indication of how strongly the sample contig
matches each reference genome.

=head2 Virtual Methods

=head3 adjust_scores

    $scoring->adjust_scores($contig);

Perform a final adjustment on the scores for the contig before forming them into a vector. Use this method
for adjustments that are related solely to the scores themselves. For adjustments that are relevant to the
vector nature of the scores, use L</adjust_vector>. The default method does nothing.

=over 4

=item contig

L<SampleContig> object whose scores are to be adjusted.

=back

=cut

sub adjust_scores {
    # Get the parameters.
    my ($self, $contig) = @_;
    # Loop through the contig's hit scores.
    my $scoreH = $contig->scores;
    for my $refID (keys %$scoreH) {
        my $hits = $contig->ref_hit($refID);
        if ($hits) {
            $scoreH->{$refID} /= $hits;
        }
    }
}

1;