package Anno;

use strict;
use Stats;
use Data::Dumper;
use SeedUtils;
use gjoseqlib;
use BlastUtils;
#######################################################

=head1 annotation_utilities

=head3 check_all_subsystems

    my @suggestions = &check_all_subsystems()

This routine checks all subsystems, returning a list
of recommended changes.

=over 4

=item RETURNS

Suggestions are returned in @suggestions and are of the form

    [PEG,change-code,optional-arguments]

In this case, they would be

    [$peg,'change-function',$function]
    ['NEW','add-peg',[$location,$function]]

and maybe a few more (e.g, deletion of a PEG)

=back

=cut

=head3 check_subsystem

    my @suggestions = &check_subsystem($subsystem)

This routine checks a single subsystem.

Suggestions are returned in @suggestions and are of the form

    [PEG,change-code,optional-arguments]

In this case, they would be

    [$peg,'change-function',$function]
    ['NEW','add-peg',[$location,$function]]

and maybe a few more.

=over 4

=item $subsystem

This argument gives access to the subsystem, from which you can
get rows in the populated subsystem.

=item RETURN

The list of suggested changes is returned.  For now, we will
simply use them to construct a report for Veronika.

=back

=cut

sub check_sbsystem {
    my($subsystem) = @_;

    my @suggestions = ();

    return @suggestions;
}

########################################

=head3 check_row_in_populated_subsystem

    my @suggestions = &check_row_in_populated_subsystem($row_id)

Checks the cells in one row of one subsystem (the subsystem and
genome ids are accessible via the row_id).  This routine checks
for empty cells and cells that contain multiple PEGs.
Suggestions are returned in @suggestions and are of the form

    [PEG,change-code,optional-arguments]

In this case, they would be

    [$peg,'change-function',$function]
    ['NEW','add-peg',[$location,$function]]

and maybe a few more.

=over 4

=item $row_id

This argument gives access to the subsystem, roles, and genome
to be checked.

=item RETURN

The list of suggested changes is returned.  For now, we will
simply use them to construct a report for Veronika.

=back

=cut

sub check_row_in_populated_sbsystem {
    my($role_id) = @_;

    my @suggestions = ();

    return @suggestions;
}

########################################

=head3 make_blastdb_for_role

    &make_blastdb($role)

Makes a blast database for the protein sequences of
all instances of $role, as well as multifunctional sequences
that implement $role.

This is redundant with respect to the blast datbases used
in building role properties.  I think that we should
have one databse per role that includes multifunctional
sequences (but, that requires that hits be filtered
when the multifunctional entries are not desired.


=over 4

=item $role

a roleid (not the description)

=item RETURN

Nothing is returned, but the constructed blastdb is saved
in the directory supporting characterization of roles.

=back

=cut

#######################################################

=head3 in_coding

    my @overlapping_pegs = &in_coding($loc)

Given a coding region (often from blast), this routine checks
to see if it is embedded in or overlaps one or more PEGs.

The code should use the L<Shrub/FeaturesInRegion> method.

=over 4

=item $loc

the location that is believed to be coding

=item RETURN

Returns a list of 3-tuples.  Each 3-tuple contains

=over 8

=item PEG id

=item PEG location

=item location of overlap

=back

=back

=cut

sub in_coding {
    my($loc) = @_;

    my @overlapping_pegs;

    return@overlapping_pegs;
}

#######################################################
#  We need
#           peg_to_subsystems($peg) [returns list of subsystems]
#
#  This should be added to Shrub.pm
#######################################################

=head3 check_length_of_role_instance

    my $z_score = check_length_of_role_instance($role,$length)

Get the z-score of the length of an instance of a role.  The length is
the length of the protein string.

=over 4

=item $role

NOTE: use the id and not the role description

=item $length

length of the protein sequence of the instance

=item RETURN

returns the z-score (the number of standard deviations from the mean -- thus,
0.0 would be rock solid)

=back

=cut

sub check_length_of_role_instance {
    my($roleid,$length) = @_;

}

#################################################


=head3 check_similarity_of_role_instance

    my $z_score = check_similarity_of_role_instance($role,$similarity)

Get the z-score of the similarity of an instance of a role.  The similarity is
the best similarity to a known instance

=over 4

=item $role

NOTE: use the id and not the role description

=item $similarity

the best normalized bit-score of a blast against known instances.  The mean
here was computed by blasting known instances against the blastdb of known instances.

=item RETURN

returns the z-score (the number of standard deviations from the mean -- thus,
0.0 would be rock solid)

=back

=cut

sub check_similarity_of_role_instance {
    my($role,$similarity) = @_;

}

#################################################

=head3 check_cdd_for_role_instance

    my $cdd = check_cdd_for_role_instance($role,$translation)

is used to see if the instance of a role (i.e., the $translation)
has a confirming cdd hit.

=over 4

=item $role

NOTE: use the id and not the role description

=item $translation

a protein string

=item RETURN

returns the z-score (the number of standard deviations from the mean -- thus,
0.0 would be rock solid)

=back

=cut

sub check_cdd_for_role_instance {
    my($role,$translation) = @_;

}

#################################################

=head3 expand_coding

    my ($expanded_loc,$truncated) = expand_coding($loc,$role,$ncbi_code)

is used to expand a coding segment to predicted start/stop codons.  We
are not attempting to be comprehensize, so we properly use only codes
11 and 4 (default is 11).  If 4 is specified TGA is not considered a STOP.

=over 4

=item $loc

The location is a region believed to be coding

=item $role (optional)

If $role is not specified, then the routine simply searches for the
START producing the longest gene.  If $role is specified, the code
searches for the start that produces a length closest to the mean
length of instances of $role.

=item RETURN

returns an attempt to find the best START and STOP codons, which may,
or may not, succeed.  If it succeeds $truncated will be FALSE (else, TRUE).
$expanded_loc gives our best guess of the location of the hypothesized peg.

=back

=cut

sub expand_coding {
    my($loc,$role,$ncbi_code) = @_;

}

#################################################

=head3 extend_to_start

    my ($expanded_loc,$truncated) = extend_to_start($loc,$role,$ncbi_code)

is used to expand a coding segment to a predicted start codon.  We
are not attempting to be comprehensize, so we properly use only codes
11 and 4 (default is 11).  If 4 is specified TGA is not considered a STOP.

=over 4

=item $loc

The location is a region believed to be coding

=item $role (optional)

If $role is not specified, then the routine simply searches for the
START producing the longest gene.  If $role is specified, the code
searches for the start that produces a length closest to the mean
length of instances of $role.

=item RETURN

returns an attempt to find the best START codon, which may,
or may not, succeed.  If it succeeds $truncated will be FALSE (else, TRUE).
$expanded_loc gives our best guess of the location of the START.

=back

=cut

sub extend_to_start {
    my($loc,$role,$ncbi_code) = @_;

}

#################################################

=head3 extend_to_stop

    my ($expanded_loc,$truncated) = extend_to_stop($loc,$ncbi_code)

is used to expand a coding segment to a predicted stop codon.  We
are not attempting to be comprehensize, so we properly use only codes
11 and 4 (default is 11).  If 4 is specified TGA is not considered a STOP.

=over 4

=item $loc

The location is a region believed to be coding

=item $ncbi_code

Only 4 and 11 are used.

=item RETURN

returns an attempt to find the best STOP codon, which may,
or may not, succeed.  If it succeeds $truncated will be FALSE (else, TRUE).
$expanded_loc gives our best guess of the location of the STOP.

=back

=cut

sub extend_to_stop {
    my($loc,$ncbi_code) = @_;

}

#################################################


=head3 broad_specificity

    my @possibilities = &broad_specificity($role)

This routine checks to see if any PEGs are believed to
implement both $role and other roles due to broad specificity.

=over 4

=item $role

A functional role ID (not a description)

=item RETURN

returns a list of possibilities.  Each possibility is a 2-tuple
containing a role ID and a PEG that is bifunctional with the
two roles.  The peg is an exemplar -- what you get back will be
a list with only one tuple for each distinct role.
The returned list may well be empty.

=back

=cut

sub broad_specificity{
    my($role) = @_;

}

#################################################

=head3 check_empty_cells

    my @suggestions = &check_empty_cells($subsystem)

Returns suggestions on filling cells that are empty in a subsystem
spreadsheet.

For each empty cell or role R and genome G, there are three main cases:

=over 4

=item 1

a blast of known instances of R finds nothing in G (no suggestions)

=item 2

a blast of known instances of R picks up a called peg that is possibly misannotated.  In this case we try to check the possibility
and make suggestions as appropriate.

=item 3

a blast of known instances picks up a probable gene that was not called. This may be really problematic,
when frameshifts exist in the potential gene.  These might be sequencing artifacts or real mutations leading to a pseudogene.

=back

In each case we attempt to unravel what should be done.

=over 4

=item $subsystem

This module handles one subsystem.  The supporting role-property data and bastdbs
should be constructed lazily.

=item RETURN

Suggestions are returned.  They amount to structured assertions of what should be done.

=back

=cut

sub check_empty_cells {
    my($subsystem) = @_;

    my @suggestions = ();

    return @suggestions;
}

#################################################

=head3 check_duplicate_cell

    my @suggestions = &check_duplicate_cell($subsystem,$role,$genome)

Returns suggestions on how to handle a cell containing multiple PEGs.
PEGs cannot be moved out of the cell if

=over 4

=item 1

it is part of a chromosomal cluster for the subsystem.

=item 2

it looks like an instance of the role in that the length and similarity values are within 3 standard deviations and
the CDD covering motiff hits within the PEG.

=back

Otherwise, pick the most likely candidate (based on z-scores of length
and similarity, along with CDD motiff coverage), and try to "move"
the remaining PEGs.

To move a PEG with function F, make the assertion that suggests
changing the assignment to

    hypothetical protein

or

    unreliable F

depending on whether or not F is hypothetical.  Note that it will
usually not be hypothetical, since it includes R, but it could
be something like "FIGxxxxx: hypothetical protein".

=over 4

=item $subsystem

The subsystem is specified to allow checking whether or not pegs are
part of a cluster (which uses roles in the subsystem and the genome
to recognize a containing cluster -- that is the clusters are not BBH
induced).

=item $role

This method handles a single cell containing multiple PEGs.  When
appropriate it makes suggestions on how to disambiguate the assignments.

=item $genome

This is any genome ID.  Genomes that do not contain multiple instances
of $role just produce no suggested actions.

=item RETURN

Suggestions are returned.  They amount to structured assertions of what should be done.

=back

=cut

sub check_duplicates_in_cell {
    my($role,$genome) = @_;

    my @suggestions = ();

    return @suggestions;
}
