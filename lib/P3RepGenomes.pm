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


package P3RepGenomes;

    use strict;
    use warnings;
    use RepGenome;
    use P3Utils;
    use RoleParse;

=head1 PATRIC Utilities for a Representative-Genome Database

This package contains services for manipulating a representative-genome database with PATRIC.

The constructor connects to both PATRIC and the representative-genome database.  The special methods use both connections
to perform representative-genome operations.

The fields in this object are as follows.

=over 4

=item repDB

The L<RepGenomeDb> object for the representative-genome database.

=item p3

The L<P3DataAPI> object for connecting to PATRIC.

=item role

The role used to search for seed proteins in the PATRIC database.

=item checksum

The checksum for the above role.

=back


=head2 Special Methods

=head3 new

    my $p3repDB = P3RepGenomes->new($p3, $repDB, %options);

Create a new PATRIC Representative Genome Database object.

=over 4

=item p3

The L<P3DataAPI> object for connecting with PATRIC.

=item repDB

The L<RepGenomeDb> object containing the representative-genome database.

=item options

A hash containing zero or more of the following options.

=over 8

=item role

The name of the role used for the seed protein.  The default is L<Phenylalanyl-tRNA synthetase alpha chain>.

=back

=back

=cut

sub new {
    my ($class, $p3, $repDB, %options) = @_;
    # Get the options.
    my $role = $options{role} // 'Phenylalanyl-tRNA synthetase alpha chain';
    # Compute the role checksum.
    my $checkSum = RoleParse::Checksum($role);
    # Create and bless the object.
    my $retVal = {
        p3 => $p3,
        repDB => $repDB,
        role => $role,
        checksum => $checkSum
    };
    bless $retVal, $class;
    # Return the object.
    return $retVal;
}


=head2 Utility Methods

=head3 FindClosest

    my $genomeHash = $p3repDB->FindClosest($prot, %options);

Find the genomes closest to an incoming genome represented by a seed protein.  The representative sets whose representative genomes are
within a specified distance to the incoming genome will be searched and the best results kept.

=over 4

=item prot

The seed protein sequence for the genome of interest.

=item options

A hash containing zero or more of the following keys.

=over 8

=item minScore

The minimum similarity score for a representative genome to be considered acceptable.  The default is C<25>.

=item N

The number of output genomes desired.  The default is C<10>.

=back

=item RETURN

Returns a reference to a hash mapping each neighboring genome ID to the pair [similarity score, genome name].

=back

=cut

sub FindClosest {
    my ($self, $prot, %options) = @_;
    # Get the RepGenome database.
    my $repDB = $self->{repDB};
    # Get the options.
    my $minScore = $options{minScore} // 25;
    my $N = $options{N} // 10;
    # Find the neighboring representative sets.
    my $repList = $repDB->list_reps($prot, $minScore);
    # Create a representative-genome object for the incoming protein.
    my $protRep = RepGenome->new('incoming', prot => $prot, K => $repDB->K);
    # Now loop through the sets, saving the closest.  This will contain the N closest genomes in the form [id, score].
    my @outList;
    # This tracks the score of the furthest genome in the result list.
    my $minFound = $minScore;
    for my $rep (keys %$repList) {
        # Search the represented genomes in this list, collecting distances to the incoming genome.
        my $repObject = $repDB->rep_object($rep);
        # Start with the representative genome itself.
        my $score = $repList->{$rep};
        _merge_genome(\@outList, \$minFound, $N, $rep, $score, $repObject->name);
        # Loop through the represented genomes.
        my @genomes = map { $_->[0] } @{$repObject->rep_list()};
        my $protHash = $self->GetProts(\@genomes);
        for my $genome (keys %$protHash) {
            my $pair = $protHash->{$genome};
            $score = $protRep->check_genome($pair->[1]);
            _merge_genome(\@outList, \$minFound, $N, $genome, $score, $pair->[0]);
        }
    }
    # Convert the results into a hash and return them.
    my %retVal = map { $_->[0] => [$_->[1], $_->[2]] } @outList;
    return \%retVal;
}

=head3 GetProts

    my $protHash = $p3repDB->GetProts(\@genomes);

Find the seed protein sequences for the specified genomes.

=over 4

=item genomes

Reference to a list of genome IDs.

=item RETURN

Returns a reference to a hash mapping each incoming genome ID to a pair [genome_name, amino acid sequence for its seed protein].

=back

=cut

sub GetProts {
    my ($self, $genomes) = @_;
    # Get the P3 object.
    my $p3 = $self->{p3};
    # Get the protein role definitions.
    my $role = $self->{role};
    my $checkSum = $self->{checksum};
    # Ask for the proteins.
    my $resultList = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', $role]], ['genome_id', 'patric_id', 'product', 'genome_name', 'aa_sequence'],
            $genomes, 'genome_id');
    # Loop through the results, keeping good ones.  Note we only keep the longest for each genome.
    my %retVal;
    for my $result (@$resultList) {
        my ($genome, $fid, $role2, $name, $seq) = @$result;
        my $checkSum2 = RoleParse::Checksum($role2);
        if ($checkSum2 eq $checkSum) {
            if (! $retVal{$genome} || length($retVal{$genome}[1]) > length($seq)) {
                $retVal{$genome} = [$name, $seq];
            }
        }
    }
    # Return the hash.
    return \%retVal;
}


=head2 Internal Utilities

=head3 _merge_genome

    _merge_genome($outList, $pMinFound, $N, $genome, $score);

Merge the specified genome and score into the output list.  The output list size is limited, and the score of the
last genome is remembered for performance enhancement.

=over 4

=item outList

A reference to the output list, which contains [genome, score, name] pairs sorted from best to worst.

=item pMinFound

A reference to a scalar containing the score of the worst genome in the list.

=item N

The desired maximum list size.

=item genome

The ID of the genome to merge.

=item score

The score of the genome to merge.

=back

=cut

sub _merge_genome {
    my ($outList, $pMinFound, $N, $genome, $score, $name) = @_;
    # Get the old list size.
    my $n = scalar @$outList;
    # Insure the new genome goes in the current list.
    if ($score <= $$pMinFound) {
        # Here we go at the end, but only if there is room.
        if ($n < $N) {
            push @$outList, [$genome, $score, $name];
            $$pMinFound = $score;
        }
    } else {
        # Here we go in the middle.
        my $i = 0;
        while ($i < $n && $outList->[$i][1] >= $score) { $i++ }
        splice @$outList, $i, 0, [$genome, $score, $name];
        # If the array is now too big, chop off the last entry.
        if ($n == $N) {
            pop @$outList;
            $$pMinFound = $outList->[$n - 1][1];
        }
    }
}

1;