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
    use FIG_Config;
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

=item uniHash

A reference to a hash keyed on universal role ID.

=item minScore

The minimum acceptable score. Lower scores are set to 0.

=item uniTotal

The total number of universal roles.

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

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item minscore

The minimum acceptable score. (Lower scores are set to 0.)

=item unifile

The name of a tab-delimited file containing the universal role IDs in the first column. The default is C<uni_roles.tbl> in
the global data directory.

=back

=back

=cut

sub new_for_script {
    my ($class, $opt) = @_;
    # Compute the universal role hash.
    my $uniHash = ReadUniHash($opt->unifile);
    # Create the object.
    my $retVal = {
        covgWeight => $opt->covgweight,
        refWeight => $opt->refweight,
        uniPenalty => $opt->unipenalty,
        uniWeight => $opt->uniweight,
        uniHash => $uniHash,
        minScore => ($opt->minscore)
    };
    # Compute the number of universal roles.
    $retVal->{uniTotal} = scalar keys %$uniHash;
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 new

    my $score = Bin::Score->new($covgweight, $refweight, $unipenalty, $uniweight, $minscore, $unifile);

Create a scoring object from specified scoring values.

=over 4

=item covgweight

The weight to assign to the coverage score.

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item minscore

The minimum acceptable score. (Lower scores are set to 0.)

=item unifile (optional)

The name of a tab-delimited file containing the universal role IDs in the first column. The default is C<uni_roles.tbl> in
the global data directory.

=back

=cut

sub new {
    my ($class, $covgweight, $refweight, $unipenalty, $uniweight, $minscore, $unifile) = @_;
    # Compute the universal role hash.
    my $uniHash = ReadUniHash($unifile);
    # Create the object.
    my $retVal = {
        covgWeight => $covgweight,
        refWeight => $refweight,
        uniPenalty => $unipenalty,
        uniWeight => $uniweight,
        uniHash => $uniHash,
        minScore => $minscore
    };
    # Compute the number of universal roles.
    $retVal->{uniTotal} = scalar keys %$uniHash;
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

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item unifile

The name of a tab-delimited file containing the universal role IDs in the first column. The default is C<uni_roles.tbl> in
the global data directory.

=item minscore

The minimum acceptable score. Lower scores are set to 0.

=back

This method returns the specifications for these command-line options in a form
that can be used in the L<ScriptUtils/Opts> method.

=cut

sub script_options {
    return (
           [ "covgweight=f",  "the weight to assign to the coverage score", { default => 0.717 } ],
           [ "refweight=f",   "the weight to assign to the closest-reference-genome score", { default => 0.255 } ],
           [ "unipenalty=f",  "the penalty to assign to duplicate universal roles", { default => 0.970 } ],
           [ "uniweight=f",   "the weight to assign to the universal role score", { default => 0.856 } ],
           [ "unifile=s",     "file containing the universal roles", { default => "$FIG_Config::global/uni_roles.tbl" } ],
           [ "minscore=f",    "the minimum acceptable score (lower scores are set to 0)", { default => 1.15 }]
    );
}

=head3 scale_min

    my $minScore = Bin::Score::scale_min($covgWeight, $refWeight, $uniWeight, $inMinScore);

Scale a minimum score in the range from 0 to 1 to be proportional to the sum of the weights. This is used to convert
a parameter chromosome from L<ga.pl> to actual scoring parameters.

=over 4

=item covgWeight

The weight to assign to the coverage score.

=item refWeight

The weight to assign to the reference genome score.

=item uniWeight

The weight to assign to the universal role score.

=item inMinScore

The unscaled minimum acceptable score.

=item RETURN

Returns the scaled minimum acceptable score.

=back

=cut

sub scale_min {
    my ($covgWeight, $refWeight, $uniWeight, $inMinScore) = @_;
    my $retVal = (1 + $inMinScore) * 0.4 * ($covgWeight + $refWeight + $uniWeight);
    return $retVal;
}

=head3 Vector

    my $vector = Bin::Score::Vector($bin1, $bin2);

Return the scoring vector for a pair of bins. The scoring vector consists of the coverage score,  the closest-reference-genome score, the number of universal roles not in common, and the number of
universal roles in common. These components are combined using the weights in the scoring object to compute the
real score.

=over 4

=item bin1

A L<Bin> object representing the first set of contigs.

=item bin2

A L<Bin> object representing the second set of contigs.

=item RETURN

Returns a 5-tuple consisting of (0) the coverage score, (1) the closest-reference-genome
score, (2) the number of universal roles not in common, and (3) the number of universal roles in common.

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


=head3 Cmp

    my $cmp = Bin::Score::Cmp($a, $b);

Compare two scoring vectors for sort purposes. The sort is by highest coverage score,
then highest reference-genome score, highest unique-universal score, lowest common-universal score.

=over 4

=item a

First scoring vector to compare.

=item b

Second scoring vector to compare.

=item RETURN

Returns a number indicating the sort order of the two vectors.

=back

=cut

sub Cmp {
    my ($a, $b) = @_;
    my $retVal = ($b->[0] <=> $a->[0]) ||
                 ($b->[1] <=> $a->[1]) ||
                 ($b->[2] <=> $a->[2]) ||
                 ($b->[3] <=> $a->[3]);
    return $retVal;
}



=head3 ReadUniHash

    my $uniHash = Bin::Score::ReadUniHash($unifile);

Read the universal role file and produce a hash keyed by universal role
ID. If no file name is provided, the default will be used.

=over 4

=item unifile

The name of a tab-delimited file containing the universal roles. The first column must contain the universal role ID and the
third should contain the role description. If this parameter is omitted, the default file C<uni_roles.tbl> in the global
data directory will be used.

=item RETURN

Returns a reference to a hash mapping each universal role ID to its description.

=back

=cut

sub ReadUniHash {
    # Get the parameters.
    my ($unifile) = @_;
    # Declare the return variable.
    my %retVal;
    # Open the file for input.
    my $file = ($unifile // "$FIG_Config::global/uni_roles.tbl");
    open(my $ih, '<', $file) || die "Could not open universal role file: $!";
    # Read the roles.
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        my ($id, undef, $name) = split /\t/, $line;
        $retVal{$id} = $name;
    }
    # Return the result.
    return \%retVal;
}


=head2 Public Manipulation Methods

=head3 Reset

    $score->Reset($covgweight, $refweight, $unipenalty, $uniweight, $minscore);

Set the scoring weights to new values.

=over 4

=item covgweight

The weight to assign to the coverage score.

=item refweight

The weight to assign to the reference genome score.

=item unipenalty

The penalty for universal roles in common.

=item uniweight

The weight to assign to the universal role score.

=item minscore

The minimum acceptable score. (Lower scores are set to 0.)

=back

=cut

sub Reset {
    my ($self, $covgweight, $refweight, $unipenalty, $uniweight, $minscore) = @_;
    $self->{covgWeight} = $covgweight;
    $self->{refWeight} = $refweight;
    $self->{uniPenalty} = $unipenalty;
    $self->{uniWeight} = $uniweight;
    $self->{minScore} = $minscore;
}

=head3 ScoreV

    my $value = $score->ScoreV(\@vector);

Compute the score from a scoring vector.

=over 4

=item vector

Vector of scores from which the final score should be computed. A 5-tuple consisting of (0) the coverage
score, (1) the closest-reference-genome score, (2) the number of universal roles not in common, and (3) the
number of universal roles in common.

=item RETURN

Returns the score computed by applying the weights to the individual scores in the vector.

=back

=cut

sub ScoreV {
    my ($self, $vector) = @_;
    my $retVal = $self->{covgWeight} * $vector->[0] + $self->{refWeight} * $vector->[1];
    my $uscore = $vector->[2] - $self->{uniPenalty} * $vector->[3];
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

=head3 uni_hash

    my $uniHash = $score->uni_hash;

Return the universal role hash.

=cut

sub uni_hash {
    my ($self) = @_;
    return $self->{uniHash};
}

=head3 uni_total

    my $uniTotal = $score->uni_total;

Return the total number of universal roles.

=cut

sub uni_total {
    my ($self) = @_;
    return $self->{uniTotal};
}

=head3 Show

    my $string = $score->Show();

Return the score parameter values as a printable string.

=cut

sub Show {
    my ($self) = @_;
    my @retVal;
    for my $parm (keys %$self) {
        my $value = $self->{$parm};
        if (! ref $value) {
            push @retVal, "$parm = $self->{$parm}";
        }
    }
    return join("\n", "** SCORING PARAMETERS", @retVal, "", "");
}


=head3 FilterBin

    $score->FilterBin($bin);

Process a bin object to remove the universal roles not active in this scoring environment.

=over 4

=item bin

A L<Bin> object to be updated.

=back

=cut

sub FilterBin {
    my ($self, $bin) = @_;
    # Get the bin's universal role hash.
    my $uniProtH = $bin->uniProts;
    # Delete the roles not in our role list.
    my $uniHash = $self->{uniHash};
    for my $role (keys %$uniProtH) {
        if (! $uniHash->{$role}) {
            delete $uniProtH->{$role};
        }
    }
}


=head3 Filter

    $score->Filter($binList);

Filter all the bins in the specified list to remove universal roles not used by this scoring method.

=over 4

=item binList

A reference to a list of L<Bin> objects to be filtered.

=back

=cut

sub Filter {
    my ($self, $binList) = @_;
    # Loop through the bins, filtering them.
    for my $bin (@$binList) {
        $self->FilterBin($bin);
    }
}


=head3 Copy

    my $newBinList = $score->Copy($binList);

Create filtered copies of all the bins in the specified list.

=over 4

=item binList

A reference to a list of L<Bin> objects to be filtered.

=item RETURN

Returns a reference to a list of filtered copies of the incoming L<Bin> objects.

=back

=cut

sub Copy {
    my ($self, $binList) = @_;
    # Create the return list by copying the bins.
    my @retVal = map { Bin->new_copy($_) } @$binList;
    # Filter the bins.
    $self->Filter(\@retVal);
    # Return the result.
    return \@retVal;
}

1;