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


package Bin::Score;

    use strict;
    use warnings;

=head1 Community Bin Scoring Object

This object is used to score the comparison between two L<Bin> objects. The six scoring parameters are stored in here,
as well as the total number of universal roles. These parameters are then used to compare two bins to determine whether
or not they belong together.

The object has the following fields.

=over 4

=item covgWeight

The weight to assign to the coverage score. The coverage score is the fraction of the coverage values that are within
20% of each other. (So, the coverage vectors are compared, and the number of coordinates in the first bin that are within
20% of the corresponding value for the second bin is computed and divided by the vector length.)

=item tetraWeight

The weight to assign to the tetranucleotide score. The tetranucleotide score is the dot product of the two tetranucleotide
vectors.

=item refWeight

The weight to assign to the reference genome score. The reference genome score is 1.0 if the two reference genome sets are
identical and nonempty, 0.6 if the two sets are empty, 0.5 if one set is empty and the other is not, and 0 if both sets
are nonempty but different.

=item uniPenalty

The penalty for universal roles in common.

=item uniWeight

The weight to assign to the universal role score. This is equal to the number of universal roles in exactly one of the two
bins less the C<uniPenalty> times the number of universal roles in both bins, all scaled by the total number of universal
roles. A negative values is changed to zero.

=item uniTotal

The total number of universal roles.

=item minScore

The minimum acceptable score. Lower scores are set to 0.

=back

=head2 Special Methods

=head3 new_for_script

    my $score = Bin::Score->new_for_script($opt);

Create a scoring object from command-line options.

=over 4

=item opt

A L<Getopt::Long::Descriptive::Opts> object containing command-line options. All of the following options must be present.

=over 8

=item covgweight

The weight to assign to the coverage score.

=item tetraweight

The weight to assign to the tetranucleotide score.

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item minscore

The minimum acceptable score. (Lower scores are set to 0.)

=item unitotal

The total number of universal roles. This is the only value with a default-- C<58>.

=back

=back

=cut

sub new_for_script {
    my ($class, $opt) = @_;
    # Create the object.
    my $retVal = {
        covgWeight => $opt->covgweight,
        tetraWeight => $opt->tetraweight,
        refWeight => $opt->refweight,
        uniPenalty => $opt->unipenalty,
        uniWeight => $opt->uniweight,
        uniTotal => ($opt->unitotal // 58),
        minScore => ($opt->minscore)
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 new

    my $score = Bin::Score->new($covgweight, $tetraweight, $refweight, $unipenalty, $uniweight, $minscore, $unitotal);

Create a scoring object from specified scoring values.

=over 4

=item covgweight

The weight to assign to the coverage score.

=item tetraweight

The weight to assign to the tetranucleotide score.

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item minscore

The minimum acceptable score. (Lower scores are set to 0.)

=item unitotal (optional)

The total number of universal roles. If omitted, C<57> is assumed.

=back

=cut

sub new {
    my ($class, $covgweight, $tetraweight, $refweight, $unipenalty, $uniweight, $minscore, $unitotal) = @_;
    # Create the object.
    my $retVal = {
        covgWeight => $covgweight,
        tetraWeight => $tetraweight,
        refWeight => $refweight,
        uniPenalty => $unipenalty,
        uniWeight => $uniweight,
        uniTotal => ($unitotal // 57),
        minScore => $minscore
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 script_options

    my @opt_specs = Bin::Score::script_options();

These are the command-line options for configuring a L<Bin::Score> object.

=over 4

=item covgweight

The weight to assign to the coverage score.

=item tetraweight

The weight to assign to the tetranucleotide score.

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item unitotal

The total number of universal roles. If omitted, C<58> is assumed.

=item minscore

The minimum acceptable score. Lower scores are set to 0.

=back

This method returns the specifications for these command-line options in a form
that can be used in the L<ScriptUtils/Opts> method.

=cut

sub script_options {
    return (
           [ "covgweight=f", "the weight to assign to the coverage score", { default => 0.25 } ],
           [ "tetraweight=f", "the weight to assign to the tetranucleotide score", { default => 0.25 }  ],
           [ "refweight=f", "the weight to assign to the closest-reference-genome score", { default => 0.25 } ],
           [ "unipenalty=f", "the penalty to assign to duplicate universal roles", { default => 1 } ],
           [ "uniweight=f", "the weight to assign to the universal role score", { default => 0.25 } ],
           [ "unitotal=i", "the total number of universal roles", { default => 57 } ],
           [ "minscore=f", "the minimum acceptable score (lower scores are set to 0)", { default => 0.5 }]
    );
}

=head3 Vector

    my $vector = Bin::Score::Vector($bin1, $bin2);

Return the scoring vector for a pair of bins. The scoring vector consists of the coverage score, the tetranucleotide
score, the closest-reference-genome score, the number of universal roles not in common, and the number of
universal roles in common. These components are combined using the weights in the scoring object to compute the
real score.

=over 4

=item bin1

A L<Bin> object representing the first set of contigs.

=item bin2

A L<Bin> object representing the second set of contigs.

=item RETURN

Returns a 5-tuple consisting of (0) the coverage score, (1) the tetranucleotide score, (2) the closest-reference-genome
score, (3) the number of universal roles not in common, and (4) the number of universal roles in common.

=back

=cut

sub Vector {
    my ($bin1, $bin2) = @_;
    my ($i, $n);
    # This will be the output vector.
    my @retVal;
    # Compare the coverage vectors.
    my $covg1 = $bin1->coverage;
    my $covg2 = $bin2->coverage;
    $n = scalar @$covg1;
    my $cScore = 0;
    for ($i = 0; $i < $n; $i++) {
        # Get the corresponding coverages.
        my ($covgV1, $covgV2) = ($covg1->[$i], $covg2->[$i]);
        # Count them if they are close.
        if ($covgV1 > $covgV2) {
            if ($covgV1 <= $covgV2 * 1.2) {
                $cScore++;
            }
        } elsif ($covgV2 > $covgV1) {
            if ($covgV2 <= $covgV1 * 1.2) {
                $cScore++;
            }
        } else {
            $cScore++;
        }
    }
    $cScore /= $n;
    push @retVal, $cScore;
    # Compare the tetranucleotide vectors.
    my $dot = 0;
    my $tetra1 = $bin1->tetra;
    my $tetra2 = $bin2->tetra;
    my $vlen = $bin1->tetraLen * $bin2->tetraLen;
    $n = scalar @$tetra1;
    for ($i = 0; $i < $n; $i++) {
        $dot += $tetra1->[$i] * $tetra2->[$i];
    }
    my $tscore = $dot / $vlen;
    push @retVal, $tscore;
    # Compare the reference genome lists.
    my @ref1 = $bin1->refGenomes;
    my @ref2 = $bin2->refGenomes;
    my $refScore;
    if (! @ref1 && ! @ref2) {
        # Both sets are empty.
        $refScore = 0.6;
    } elsif (! @ref1 || ! @ref2) {
        # One set is empty.
        $refScore = 0.5;
    } else {
        # Both sets are nonempty. Compare the sets.
        $n = scalar @ref1;
        if ($n != scalar(@ref2)) {
            $refScore = 0;
        } else {
            # Sets are the same size. The elements are sorted, so we can
            # do a straight compare.
            $refScore = 1;
            for ($i = 0; $i < $n && $refScore; $i++) {
                if ($ref1[$i] ne $ref2[$i]) {
                    $refScore = 0;
                }
            }
        }
    }
    push @retVal, $refScore;
    # Now check the universal proteins. We track the number in both bins and the number in only one bin.
    my $uCommon = 0;
    my $uOnly = 0;
    my $univ1 = $bin1->uniProts;
    my $univ2 = $bin2->uniProts;
    my $u;
    for $u (keys %$univ1) {
        if ($univ2->{$u}) {
            $uCommon++;
        } else {
            $uOnly++;
        }
    }
    for $u (keys %$univ2) {
        if (! $univ1->{$u}) {
            $uOnly++;
        }
    }
    push @retVal, $uOnly, $uCommon;
    # Return the scoring vector.
    return \@retVal;
}


=head2 Public Manipulation Methods

=head3 ScoreV

    my $value = $score->ScoreV(\@vector);

Compute the score from a scoring vector.

=over 4

=item vector

Vector of scores from which the final score should be computed. A 5-tuple consisting of (0) the coverage
score, (1) the tetranucleotide score, (2) the closest-reference-genome score, (3) the number of universal
roles not in common, and (4) the number of universal roles in common.

=item RETURN

Returns the score computed by applying the weights to the individual scores in the vector.

=back

=cut

sub ScoreV {
    my ($self, $vector) = @_;
    my $retVal = $self->{covgWeight} * $vector->[0] + $self->{tetraWeight} * $vector->[1] + $self->{refWeight} * $vector->[2];
    my $uscore = $vector->[3] - $self->{uniPenalty} * $vector->[4];
    if ($uscore < 0) {
        $uscore = 0;
    } else {
        $uscore /= $self->{uniTotal};
    }
    $retVal += $self->{uniWeight} * $uscore;
    if ($retVal < $self->{minScore}) {
        $retVal = 0;
    }
    return $retVal;
}


=head3 Score

    my $value = $score->Score($bin1, $bin2);

Compute the comparison score for two bins.

=over 4

=item bin1

A L<Bin> object representing the first set of contigs.

=item bin2

A L<Bin> object representing the second set of contigs.

=item RETURN

Returns a score comparing the two bins. A high score means the bins are more likely to belong together.

=back

=cut

sub Score {
    my ($self, $bin1, $bin2) = @_;
    my $vector = Vector($bin1, $bin2);
    return $self->ScoreV($vector);
}


=head3 Show

    my $string = $score->Show();

Return the score parameter values as a printable string.

=cut

sub Show {
    my ($self) = @_;
    my @retVal;
    for my $parm (keys %$self) {
        push @retVal, "$parm = $self->{$parm}";
    }
    return join("\n", "** SCORING PARAMETERS", @retVal, "", "");
}

1;