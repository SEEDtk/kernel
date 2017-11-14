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


package RolesInSubsystems;

    use strict;
    use warnings;
    use RoleParse;
    use Carp;

=head1 Manage Roles In Subsystems

This object manages roles and role IDs. It is initialized from a roles.in.subsystems file, which is a tab-delimited file containing
(0) a role ID, (1) a role checksum, (2) a role name, and (3) an optional C<X> indicating an experimental role.
Additional roles can be added and these are marked as non-subsystem. The fields in this object are as follows.

=over 4

=item checkMap

Reference to a hash mapping checksums to role IDs.

=item nameMap

Reference to a hash mapping role IDs to role names.

=item expMap

Reference to a hash containing the IDs of experimental roles.

=item subMap

Reference to a hash containing the IDs of subsystem roles (this includes experimental roles).

=item shrub

A L<Shrub> object for accessing the database.

=back

=head2 Special Methods

=head3 new

    my $roleMap = RolesInSubsystems->new($shrub, $fileName);

Create a new RolesInSubsystems object populated with subsystem roles.

=over 4

=item shrub

A L<Shrub> object for accessing the database.

=item fileName

The name of a roles.in.subsystems file. This is a tab-delimited file, each record containing (0) a role ID, (1) a role checksum, (2) a role name,
and (3) an option C<X>, indicating an experimental role.

=back

=cut

sub new {
    my ($class, $shrub, $fileName) = @_;
    # Create the maps.
    my (%checkMap, %nameMap, %expMap, %subMap);
    # Create the object.
    my $retVal = {
        shrub => $shrub,
        checkMap => \%checkMap,
        nameMap => \%nameMap,
        expMap => \%expMap,
        subMap => \%subMap
    };
    # Read the file.
    open(my $ih, '<', $fileName) || die "Could not open $fileName: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        my ($id, $checksum, $name, $xFlag) = split /\t/, $line;
        $checkMap{$checksum} = $id;
        $nameMap{$id} = $name;
        $subMap{$id} = 1;
        if ($xFlag) {
            $expMap{$id} = 1;
        }
    }
    # Bless and return the object.
    bless $retVal, $class;
    return $retVal;
}

=head2 Public Manipulation Methods

=head3 RoleID

    my $id = RoleID($name);

Determine the ID for a role from the role name.

=over 4

=item name

Role description string for the role in question.

=item RETURN

Returns the role ID, or C<undef> if the role is not in the database.

=back

=cut

sub RoleID {
    my ($self, $name) = @_;
    # Get the checksum hash.
    my $checkMap = $self->{checkMap};
    # Compute the role checksum.
    my $checksum = RoleParse::Checksum($name);
    # Try to find the ID.
    my $retVal = $checkMap->{$checksum};
    if (! $retVal) {
        # Here we need to do a lookup.
        my $shrub = $self->{shrub};
        ($retVal) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checksum], 'id');
        if ($retVal) {
            # We found the role, so cache it for the future.
            $checkMap->{$checksum} = $retVal;
            $self->{nameMap}{$retVal} = $name;
        }
    }
    # Return the ID found.
    return $retVal;
}

=head2 Query Methods

=head3 in_sub

    my $subFlag = $roleMap->in_sub($roleID);

Return TRUE if the role is in a subsystem, else FALSE.

=over 4

=item roleID

ID of the role in question.

=item RETURN

Returns TRUE if the role is in a subsystem, else FALSE.

=back

=cut

sub in_sub {
    my ($self, $roleID) = @_;
    my $retVal = $self->{subMap}{$roleID};
    return $retVal;
}


=head3 is_exp

    my $expFlag = $roleMap->is_exp($roleID);

Return C<X> if the role is experimental, else an empty string.

=over 4

=item roleID

ID of the role in question.

=item RETURN

Returns C<X> if the role is experimental, and an empty string if it is not.

=back

=cut

sub is_exp {
    my ($self, $roleID) = @_;
    my $retVal = ($self->{expMap}{$roleID} ? 'X' : '');
    return $retVal;
}

=head3 name

    my $name = $roleMap->name($roleID);

Return the name of the specified role. It is an error if the role name is unknown.

=over 4

=item roleID

ID of the role in question.

=item RETURN

Returns the description string of the identified role.

=back

=cut

sub name {
    my ($self, $roleID) = @_;
    my $retVal = $self->{nameMap}{$roleID};
    if (! $retVal) {
        confess "Role $roleID not found in role map.";
    }
    return $retVal;
}

1;