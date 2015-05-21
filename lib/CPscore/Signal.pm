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


package CPscore::Signal;

    use strict;
    use warnings;
    use base qw(CPscore);

=head1 Community Pipeline Vector Scoring

This is a community pipeline scoring object that computes the similarity between two sample contigs based on
similarity signal vectors. In this case, the similarity vector contains the sum of the percent identity times the length
of all the BLAST matches to each reference genome. The resulting number is an indication of how much of the sample contig
matches each reference genome.

=head2 Virtual Methods

=head3 compute_score

    my $score = $scoring->compute_score($contig, $sim, $refID, $oldScore);

Process a similarity and use it to update the scoring vector coordinates for the relevant genome.

=over 4

=item contig

A L<SampleContig> object containing the data for the contig.

=item sim

A L<Sim> object containing the results of a significant BLAST hit between the contig and a
reference genome.

=item refID

The ID of the target reference genome.

=item oldScore

The current score for the similarity between the contig and the reference genome.

=back

=cut

sub compute_score {
    # Get the parameters.
    my ($self, $contig, $sim, $refID, $oldScore) = @_;
    # Get the match strength.
    my $strength = $sim->iden * abs($sim->e2() - $sim->b2()) / 100;
    my $retVal = $oldScore + $strength;
    return $retVal;
}


1;