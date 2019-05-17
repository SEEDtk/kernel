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


package TaxCheck;

    use strict;
    use warnings;
    use FIG_Config;
    use P3Utils;
    use Bin::Blast;
    use Loader;


=head1 Estimate Contig Taxonomy

This object allows processing of contig files to estimate taxonomic information for sets of contigs.  It is styled to permit
mass operation:  the tuning parameters are specified in a constructor, and the files are processed individually.

The script locates a PheS protein in the genome, then matches it to the closest good genome in PATRIC and infers the
appropriate taxonomic grouping.  Experience has shown this process is very accurate at the genus level, but less so at
more detailed levels.

The fields in this object are as follows.

=over 4

=item debug

If TRUE, progress messages will be written to STDERR.

=item protFile

A FASTA file containing examples of the universal role to use for seeding the bin assignment.

=item seedFastaFile

The name of the BLAST database for the seed protein in the various PATRIC genomes.

=item maxE

The maximum acceptable E-value.

=item refMaxE

The maximum acceptable E-value for blasting to determine the best reference genome for a seed contig. Each seed
contig eventually forms a bin.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>).

=item p3

The L<p3DataAPI> object used to access PATRIC.

=item taxCache

A reference to a hash mapping each genome ID from the main seed protein database to a taxonomic lineage.  The
lineage will be a list of 3-tuples, from largest rank to smallest, each tuple consisting of (0) a taxon ID, (1)
a taxon rank, and (2) the taxon name.

=back

=head2 Special Methods

=head3 new

    my $taxChecker = TaxCheck->new($p3, %options);

Create a new taxonomy-checker object that can be used for taxonomy estimation on contig FASTA files.

=over 4

=item p3

A L<P3DataAPI> object for accessing the PATRIC database.

=item options

A hash containing zero or more of the following keys.

=over 8

=item debug

If TRUE, progress messages will be written to STDERR.

=item protFile

A FASTA file containing examples of the universal role to use for seeding the bin assignment.  The default is
C<seedprot.fa> in the global data directory.

=item seedFastaFile

The name of the BLAST database for the seed protein in the various PATRIC genomes. The default is
C<PhenTrnaSyntAlph.fa> in the global data directory.

=item maxE

The maximum acceptable E-value. The default is C<1e-20>.

=item refMaxE

The maximum acceptable E-value for blasting to determine the best reference genome for a seed contig. Each seed
contig eventually forms a bin. The default is C<1e-10>.

=item gap

The maximum permissible gap between BLAST hits that are to be merged. BLAST hits on the same contig in the same
direction that are closer than this number of base pairs are merged into a single hit. The default is C<600>.

=item minlen

The minimum fraction length for a BLAST hit. A BLAST hit that matches less than this fraction of a protein's
length will be discarded. This is done after the gap-merging (see C<gap>). The default is C<0.50>.

=back

=back

=cut

sub new {
    my ($class, $p3, %options) = @_;
    # Extract the options.
    my $debug = $options{debug} // 0;
    my $protFile = $options{protFile} // "$FIG_Config::p3data/seedprot.fa";
    my $seedFastaFile = $options{seedFastaFile} // "$FIG_Config::p3data/PhenTrnaSyntAlph.fa";
    my $maxE = $options{maxE} // 1e-20;
    my $refMaxE = $options{refMaxE} // 1e-10;
    my $gap = $options{gap} // 600;
    my $minlen = $options{minlen} // 0.50;
    # Create the loader.
    my $loader = Loader->new();
    # Form the object.
    my $retVal = {
        debug => $debug,
        protFile => $protFile,
        seedFastaFile => $seedFastaFile,
        maxE => $maxE,
        refMaxE => $refMaxE,
        gap => $gap,
        minlen => $minlen,
        p3 => $p3,
        taxCache => {},
        loader => $loader,
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Public Methods

=head3 Compute

    my ($id, $name, \@lineage) = $taxChecker->Compute($fastaFile, $rank);

Estimate the taxonomic grouping for a FASTA file at the specified lineage ranking.

=over 4

=item fastaFile

The contig FASTA file representing a genome whose taxonomic grouping needs to be estimated.

=item rank

The desired level of granularity, either C<species> or C<genus>.

=item RETURN

Returns a three-element list containing (0) the ID of the desired taxonomic grouping, (1) the name of the desired
taxonomic grouping, and (2) a reference to a list containing the names of the groups in the full lineage, in order
from largest to smallest.

=back

=cut

sub Compute {
    my ($self, $fastaFile, $rank) = @_;
    my $maxE = $self->{maxE};
    my $minlen = $self->{minlen};
    my $gap = $self->{gap};
    my $protFile = $self->{protFile};
    my $debug = $self->{debug};
    my $refMaxE = $self->{refMaxE};
    my $seedFastaFile = $self->{seedFastaFile};
    my $loader = $self->{loader};
    # These will be the return values.
    my ($retID, $retName);
    my $retLineage = [];
    # Create the blaster.
    my $blaster = Bin::Blast->new($FIG_Config::temp, $fastaFile,
            maxE => $maxE, minlen => $minlen, gap => $gap, silent => 1);
    # Find the seed protein.
    print STDERR "Searching for seed protein in $fastaFile using $protFile.\n" if $debug;
    my $matches = $blaster->FindProtein($protFile);
    # Only proceed if we found a single match.
    my ($contig, @others) = keys %$matches;
    if (! $contig) {
        print STDERR "No seed proteins found.\n" if $debug;
    } elsif (scalar @others) {
        print STDERR "Too many seed proteins found.\n" if $debug;
    } else {
        # Get the DNA for the protein.
        my $seqHash = $loader->GetDNA($matches, $fastaFile);
        # Get the best match for the DNA.
        print STDERR "Searching for best DNA match to seed protein in $contig using $seedFastaFile.\n" if $debug;
        my $contigHash = $blaster->MatchProteins($seqHash, undef, 1, $refMaxE, db => $seedFastaFile, type => 'dna');
        # MatchProteins is a general-purpose thing that returns a list of genomes for each contig. The "1" in the parameter
        # list insures that it returns at most one.
        my $genomeData = $contigHash->{$contig}[0];
        if (! $genomeData) {
            print STDERR "No matching genome found for $contig.\n" if $debug;
        } else {
            my ($genome, $score, $name) = @$genomeData;
            # Compute the lineage for this genome.
            my $lineageL = $self->_FindLineage($genome);
            # Only proceed if we found something.
            if ($lineageL) {
                # Search for the desired rank.
                for my $tuple (@$lineageL) {
                    my ($id1, $rank1, $name1) = @$tuple;
                    # Save the taxon name.
                    push @$retLineage, $name1;
                    if ($rank1 eq $rank) {
                        # Save this ID and name if it is the appropriate rank.
                        ($retID, $retName) = ($id1, $name1);
                    }
                }
            }
        }
    }
    # Return the results found.
    return ($retID, $retName, $retLineage);
}


=head2 Internal Methods

=head3 _FindLineage

    my \@lineage = $taxChecker->_FindLineage($genome);

Return the lineage list for the specified genome ID.  The lineage will be cached in the internal hash in case it is needed again.

=over 4

=item genome

The ID of the genome whose lineage is desired.

=item RETURN

Returns a reference to a list of 3-tuples for the taxonomic groupings in the lineage, from biggest to smallest.  Each 3-tuple
consists of (0) a taxonomic group ID, (1) a rank, and (2) a taxonomic group name.

=back

=cut

sub _FindLineage {
    my ($self, $genome) = @_;
    my $p3 = $self->{p3};
    my $debug = $self->{debug};
    my $taxCache = $self->{taxCache};
    # Check the cache first.
    my $retVal = $taxCache->{$genome};
    if (! $retVal) {
        # Get the genome's lineage.
        print STDERR "Retrieving taxonomic lineage from $genome.\n" if $debug;
        my $genomeRecords = P3Utils::get_data_keyed($p3, genome => [], ['taxon_lineage_ids'], [$genome]);
        if (! @$genomeRecords) {
            print STDERR "Genome $genome not found in PATRIC.\n" if $debug;
        } else {
            my $taxonList = $genomeRecords->[0][0];
            if (! $taxonList) {
                print STDERR "Lineage not available for $genome.\n" if $debug;
            } else {
                # Get the name and rank for each element of the lineage.
                print STDERR "Retrieving lineage data.\n" if $debug;
                my $taxonRecords = P3Utils::get_data_keyed($p3, taxonomy => [],
                    ['taxon_id', 'taxon_rank', 'taxon_name'], $taxonList);
                # Pop off the over-domain stuff.
                my $n = scalar(@$taxonRecords) - 1;
                while ($taxonRecords->[$n][2] =~ /^(?:root|cellular)/i) {
                    pop @$taxonRecords;
                    $n--;
                }
                # Save the rest in a hash.
                my %taxonHash = map { $_->[0] => $_ } @$taxonRecords;
                # Roll the lineage into the output list.
                for my $taxID (@$taxonList) {
                    my $tuple = $taxonHash{$taxID};
                    if ($tuple) {
                        push @$retVal, $tuple;
                    }
                }
            }
        }
    }
    return $retVal;
}


1;