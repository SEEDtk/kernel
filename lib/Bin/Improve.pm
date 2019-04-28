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


package Bin::Improve;

    use strict;
    use warnings;
    use GEO;
    use FastA;
    use EvalCon;
    use Stats;
    use P3DataAPI;

=head1 Use Quality Information to Produce an Improved Bin

This object operates on a L<Bin> object for which a FASTA file and a C<GenomeTypeObject> exists.  The latter must have come from
RASTtk, so that it contains quality information.  We will search for contigs with no good roles and remove them, creating a
new FASTA file.

This object contains the following fields.

=over 4

=item workDir

The name of the working directory containing the bins.

=item roleHashes

A 2-tuple containing (0) a map from each role ID to its name and (1) a map from each role checksum to a role ID.

=item p3

A L<P3DataAPI> object for accessing PATRIC.

=item stats

A L<stats> object for tracking statistics.

=item minComplete

The minimum completeness for a genome to be eligible for improvement.

=back

=head2 Special Methods

=head3 new

    my $improver = Bin::Improve->new($workDir, %options);

Create a new bin improvement object.

=over 4

=item workDir

The name of the working directory containing the binning files.

=item options

A hash containing zero or more of the following options.

=over 8

=item p3

A L<P3DataAPI> object for accessing PATRIC.  If none is specified, one will be created.

=item stats

A L<Stats> object for tracking statistics.  If none is specified, one will be created.

=item minComplete

The minimum completeness for a genome to be eligible for improvement.  The default is C<90>.

=back

=back

=cut

sub new {
    my ($class, $workDir, %options) = @_;
    # Get the helper objects.
    my $p3 = $options{p3} // P3DataAPI->new();
    my $stats = $options{stats} // Stats->new();
    # Get the tuning options.
    my $min = $options{minComplete} // 90;
    # Get the role hashes.
    my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::global/roles.in.subsystems", $stats);
    # Create the object.
    my $retVal = { workDir => $workDir, roleHashes => [$nMap, $cMap], stats => $stats, p3 => $p3, minComplete => $min };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Query Methods

=head3 eligible

    my $okFlag = $improver->eligible($gto);

Return TRUE if the specified genome is eligible for improvement, else FALSE.

=over 4

=item gto

A L<GenomeTypeObject> representing the genome to be checked.  It must include quality information.  (This will be the case if it
is output from RASTtk.)

=item RETURN

Return TRUE if the genome is eligible, else FALSE.

=back

=cut

sub eligible {
    my ($self, $gto) = @_;
    my $retVal;
    my $stats = $self->{stats};
    my $quality = $gto->{quality};
    if (! $quality) {
        $stats->Add(improveNoData => 1);
    } elsif ($quality->{completeness} >= $self->{minComplete} && ! GEO::contamX($quality->{contamination})) {
        $retVal = 1;
    } else {
        $stats->Add(improveNotEligible => 1);
    }
    return $retVal;
}


=head2 Public Manipulation Methods

=head3 Improve

    my $triples = $improver->Improve($bin, $gto, $fastaFile);

Attempt to improve a bin.

=over 4

=item bin

A L<Bin> object for the bin to be improved.

=item gto

The L<GenomeTypeObject> containing the genome and its quality information.

=item fastaFile

The name to give to the FASTA file for the bin.

=item triples

Reference to a list of 3-tuples representing the contigs in the genome.

=item RETURN

Returns a reference to a list of FASTA triples for the improved bin, or C<undef> if the bin is not improvable.

=back

=cut

sub Process {
    my ($self, $bin, $gto, $fastaFile, $triples) = @_;
    # Get the stats object.
    my $stats = $self->{stats};
    # Get the work directory.
    my $workDir = $self->{workDir};
    # This will be set to a list of output triples if we improve the bin.
    my $retVal;
    # Create the GEO options.
    my %gOptions = (roleHashes => $self->{roleHashes}, detail => 2, p3 => $self->{p3}, stats => $stats);
    # Create the GEO for the sample bin.
    my $geo = GEO->CreateFromGto($gto, %gOptions);
    # Get the reference genome IDs.
    my @refGenomes = $bin->refGenomes;
    for my $refGenome (@refGenomes) {
        my $refGeo;
        if (-s "$workDir/$refGenome.json") {
            $refGeo = GEO->CreateFromGto("$workDir/$refGenome.json", %gOptions);
        } else {
            my $gHash = GEO->CreateFromPatric($refGenome, %gOptions);
            $refGeo = $gHash->{$refGenome};
        }
        $geo->AddRefGenome($refGeo);
    }
    # Do a full quality analysis.
    $geo->AnalyzeQualityData();
    # Find the bad contigs.
    my $badHash = $geo->FindBadContigs();
    my $badFound = scalar keys %$badHash;
    if (! $badFound) {
        $stats->Add(improveNoBadContigs => 1);
    } else {
        # Write and save the good contigs.
        $retVal = [];
        open(my $oh, '>', $fastaFile) || die "Could not open FASTA file $fastaFile: $!";
        for my $contig (@$triples) {
            if ($badHash->{$contig->[0]}) {
                $stats->Add(improveBadContig => 1);
            } else {
                print $oh ">$contig->[0]\n$contig->[2]\n";
                push @$retVal, $contig;
                $stats->Add(improveGoodContig => 1);
            }
        }
    }
    return $retVal;
}


1;


