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


package CPscore::Basis::HotGroup;

    use strict;
    use warnings;
    use base qw(CPscore::Basis);

=head1 Basis Vector Computation -- Strongest Signals Method

This object is used to compute the basis vector for a community sample scoring operation. The basis vector
is the list of reference genomes whose coordinates are to be used. In this case, a complicated algorithm is
used to find the reference genomes with the strongest signals.  We start by sorting each contig's scores
to find the reference genomes with the highest hits. (By default, we will look at the top four.) We then
iterate over the sorted lists, choosing the reference genome with the most top-hit positions, then
removing all the contigs in which the genome occurs at any position. The effect is to find the genome with
the strongest impact, then remove any contigs which it impacts and search for the genome with the strongest
impact on the remaining contigs, and so forth.

In addition to the fields in the base class object, this object contains the following.

=over 4

=item topSize

The number of reference genomes considered to be indicative of high impact. Only this number of genomes
will be retained when forming the list of strongest signals for a contig. The default is C<4>.

=back

=head2 Special Methods

=head3 new

    my $basis = CPscore::Basis::HotGroup->new(%options);

Construct a new basis-computation object.

=over 4

=item options

A hash of options with the following keys.

=over 8

=item topSize

The number of reference genomes considered to be indicative of high impact. Only this number of genomes
will be retained when forming the list of strongest signals for a contig. The default is C<4>.

=back

=cut

sub new {
    # Get the parameters.
    my ($class, %options) = @_;
    # Compute the options.
    my $topSize = $options{topSize} // 4;
    # Create the object.
    my $retVal = { topSize => $topSize };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head2 Virtual Methods

=head3 type

    my $typeString = $basis->type;

Return a string description of this basis computation method.

=cut

sub type {
    my ($self) = @_;
    return "hot-$self->{topSize}";
}

=head3 compute

    my $genomeList = $basis->compute(\%contigHash, \@genomeIDs);

Compute the list of reference genomes to use as the basis for the scoring vectors. This will be a subset
of the incoming genome ID list.

=over 4

=item contigHash

Reference to a hash of L<SampleContig> objects containing the scoring data for all the sample contigs,
keyed by contig ID.

=item genomeIDs

Reference to a list of genome IDs. Only genome IDs in this list can be part of the basis.

=item RETURN

Returns a list of genome IDs whose scores should be represented in the scoring vectors.

=back

=cut

sub compute {
    # Get the parameters.
    my ($self, $contigHash, $genomeIDs) = @_;
    # Create a list of contig vectors. Each vector consists of the IDs of the N
    # best-hitting genomes (in order), where N is the topSize parameter.
    my $lastI = $self->{topSize} - 1;
    my @contigList;
    # This hash tracks the number of times each genome was the best.
    my %hits;
    # Loop through the contigs.
    for my $contigID (keys %$contigHash) {
        # Get the scoring hash for this contig.
        my $contig = $contigHash->{$contigID};
        my $scores = $contig->scores;
        # Sort out the best genomes.
        my @sorted = sort { $scores->{$b} <=> $scores->{$a} } keys %$scores;
        # If there were any scores, store them. Note the actual contig ID is irrelevant.
        if (@sorted) {
            if ($#sorted < $lastI) {
                push @contigList, \@sorted;
            } else {
                push @contigList, [ @sorted[0 .. $lastI ] ];
            }
            $hits{$sorted[0]}++;
        }
    }
    # This will be our return vector.
    my @retVal;
    # Loop until we've run out of contigs to check.
    while (scalar(@contigList) > 0) {
        # Compute the winning genome.
        my ($winner, $winCount) = (undef, 0);
        for my $genomeID (keys %hits) {
            if ($hits{$genomeID} > $winCount) {
                $winner = $genomeID;
                $winCount = $hits{$genomeID};
            }
        }
        push @retVal, $winner;
        if (! defined $winner) {
            print STDERR "Winner failure.\n"; ##TODO remove this section
        }
        # When we're done, this list will contain all the contigs that didn't contain the
        # last winning genome.
        my @newContigList;
        for my $contigData (@contigList) {
            if (grep { $_ eq $winner } @$contigData) {
                # Here the winner is in this contig, so remove its influence from the
                # hit counts.
                $hits{$contigData->[0]}--;
            } else {
                # Here the winner is not in this contig, so push it into the keep list.
                push @newContigList, $contigData;
            }
        }
        # Save the new contig list and repeat.
        @contigList = @newContigList;
    }
    # Return the computed list of genome IDs;
    return \@retVal;
}

1;