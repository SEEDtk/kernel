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


package GeoGroup;

    use strict;
    use warnings;
    use GEO;
    use EvalCon;
    use Stats;
    use P3DataAPI;


=head1 GEO Genome Group

This object represents a set of L<GEO> objects that define a group of genomes.  It contains useful utilities for
managing the options and loading from directories.

The fields in this object are as follows.

=over 4

=item optionsH

Reference to a hash of the options to use in constructing the GEOs.

=item gHash

Reference to a hash mapping genome IDs to GEO objects.

=item logH

Reference to an open file handle for status messages, or C<undef> if no status messages are desired.

=item stats

Reference to a L<Stats> object for statistics.

=back

=head2 Special Methods

=head3 new

    my $geoGroup = GeoGroup->new(\%options, $gtoDir);

Create a new, empty GEO group.

=over 4

=item options

A reference to a hash of options to use in constructing the L<GEO> objects.  Missing options (C<stats>, C<p3>,
C<roleHashes>) will be filled in automatically from defaults.  This will modify the hash in place, so that it can be
used in constructors for related groups.

=item gtoDir (optional)

The name of a directory containing L<GenomeTypeObject> files to load into the group.  If omitted, the group will
remain empty.

=back

=cut

sub new {
    my ($class, $options, $gtoDir) = @_;
    # Fill in the missing options.
    my $stats = $options->{stats} // Stats->new();
    $options->{stats} //= $stats;
    my $logH = $options->{logH};
    my $p3 = $options->{p3} // P3DataAPI->new();
    $options->{p3} //= $p3;
    if (! $options->{roleHashes}) {
        my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::global/roles.in.subsystems", $stats);
        $options->{roleHashes} = [$nMap, $cMap];
    }
    # Create the object.
    my $retVal = { optionsH => $options, gHash => {}, logH => $logH, stats => $stats};
    # Bless the object.
    bless $retVal, $class;
    # Do we have a GTO directory?
    if ($gtoDir) {
        # Yes, load it.
        $retVal->LoadGtoDirectory($gtoDir);
    }
    return $retVal;
}

=head2 Public Manipulation Methods

=head3 LoadGtoDirectory

    my $count = $geoGroup->LoadGtoDirectory($gtoDir);

Load all the genomes with GTO files in the specified directory into this object.

=over 4

=item gtoDir

The name of the directory containing the GTO files to load.

=item RETURN

Returns the number of GTOs read.

=back

=cut

sub LoadGtoDirectory {
    my ($self, $gtoDir) = @_;
    # Get the stats object and the debug handle.
    my $debug = $self->{logH};
    my $stats = $self->{stats};
    # Get the input directory and extract the GTO file names.
    opendir(my $dh, $gtoDir) || die "Could not open GTO directory $gtoDir: $!";
    my @files = map { "$gtoDir/$_" } grep { $_ =~ /\.gto$/ } readdir $dh;
    closedir $dh;
    $stats->Add(gtoDir => 1);
    # Read the GTOs into GEOs.
    my $options = $self->{optionsH};
    my $newGHash = GEO->CreateFromGtoFiles(\@files, %$options);
    my $oldGHash = $self->{gHash};
    # Copy them into the hash.
    my $retVal = 0;
    for my $newGenome (keys %$newGHash) {
        $oldGHash->{$newGenome} = $newGHash->{$newGenome};
        $stats->Add(gtoAdded => 1);
        $retVal++;
    }
    # Return the count;
    return $retVal;
}

=head2 Query Methods

=head3 MapGroups

    my (\@pairs, \@orphans1, \@orphans2) = $geoGroup->MapGroups($group2, $min);

Map the genomes in this group with the genomes in another group using bidirectional best hits for the seed
proteins.

=over 4

=item group2

A L<GeoGroup> object containing the group to map to this one.

=item min (optional)

The minimum acceptable score for a match.  If omitted, it defaults to C<0>.

=item RETURN

Returns a three-element list, consisting of (0) a reference to a list of 3-tuples describing the mapped pairs in the
form [id1, id2, score], (1) a reference to a list of the IDs in the this group not having a match, and (2) a reference
to a list of the IDs in the second group not having a match.

=back

=cut

sub MapGroups {
    my ($self, $group2, $min) = @_;
    # We need to build the seed tables for the two groups.
    my $prots = $self->_seedHash;
    my $prots2 = $group2->_seedHash;
    # Compute the mappings.
    my ($pairs, $orphans, $orphans2) = GEO::FindBBHs($prots, $prots2, $min);
    # Return the results.
    return ($pairs, $orphans, $orphans2);
}

=head3 options

    my $optionsH = $geoGroup->options;

Return the options hash in this object for creating new L<GEO> objects.

=cut

sub options {
    my ($self) = @_;
    return $self->{optionsH};
}

=head3 stats

    my $stats = $geoGroup->stats;

Return the statistics object.

=cut

sub stats {
    my ($self) = @_;
    return $self->{stats};
}

=head3 geo

    my $geo = $geoGroup->geo($genome);

Return the L<GEO> for the specified genome ID.

=over 4

=item genome

The ID of the genome desired.

=item RETURN

Returns a L<GEO> for the desired genome, or C<undef> if the genome is not in this group.

=back

=cut

sub geo {
    my ($self, $genome) = @_;
    my $retVal = $self->{gHash}{$genome};
    return $retVal;
}

=head3 genomes

    my $genomesL = $geoGroup->genomes;

Return a reference to a list of the IDs for the genomes in this group.

=cut

sub genomes {
    my ($self) = @_;
    my @retVal = sort keys %{$self->{gHash}};
    return \@retVal;
}

=head3 geoList

    my $geoL = $geoGroup->geoList;

Return a reference to a list of the L<GEO> objects in this group.

=cut

sub geoList {
    my ($self) = @_;
    my @retVal;
    my $gHash = $self->{gHash};
    for my $genome (sort keys %$gHash) {
        push @retVal, $gHash->{$genome};
    }
    return \@retVal;
}

=head2 Internal Utilities

=head3 _seedHash

     my \%prots = $geoGroup->_seedHash();

This method returns a reference to a hash mapping all the genome IDs for the genomes in this group to their seed
protein sequences.

=cut

sub _seedHash {
    my ($self) = @_;
    # This will be the return value.
    my %retVal;
    # Loop through the genomes, collecting seed proteins.
    my $gHash = $self->{gHash};
    for my $genome (keys %$gHash) {
        $retVal{$genome} = $gHash->{$genome}->seed;
    }
    return \%retVal;
}

1;


