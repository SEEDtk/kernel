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


package Signatures::Raw;

    use strict;
    use warnings;
    use base qw(Signatures);

=head1 Signatures Class with Scaled Scoring

This is a subclass of L<Signatures> that uses the scaled scoring algorithm. The total score is

    hits * 1000 / (positions * groupCount)

=head2 Virtual Overrides

=head3 score

    my $score = $sigsObject->score($hits, $groupCount, $positions);

Compute the score for this group given the specified number of hits.

=over 4

=back

=cut

sub score {
    my ($self, $hits, $groupCount, $positions) = @_;
    my $retVal = $hits * 1000 / ($positions * $groupCount);
    return $retVal;
}

1;