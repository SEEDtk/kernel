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


package Bin::Blast;

    use strict;
    use warnings;
    use gjo::BlastInterface;
    use SeedUtils;
    use Stats;

=head1 Blast Analysis Object

This object BLASTs a contig against the proteins in a set of reference genomes. It returns the two closest
reference genomes and the list of universal roles found. A reference genome is only considered close if there
is a successful BLAST of a certain minimum length (default 160) and percent identity (default 0.9). The genome
with the highest percent identity match is considered the closest.

During initialization, this object builds a BLAST database of the reference genomes in a specified working
directory. If that database already exists, the initialization is significantly faster.

The fields in this object are as follows.

=over 4

=item workDir

Name of the working directory.

=item minsim

The minimum acceptable similarity for a BLAST result.

=item minlen

The minimum acceptable length for a BLAST result.

=item stats

Statistics about the blast database. This is a hash with the following values.

=over 8

=item genomes

Number of reference genomes.

=item fids

Number of proteins.

=item uniRoles

Number of universal roles.

=item uniFids

Number of proteins with universal roles.

=back

=back

In the working directory, there will be a subdirectory created for each reference genome. In that subdirectory,
The FASTA file containing the proteins will be C<blast_ref.fa>. There will also be a universal role hash in JSON format
called C<blast_ref.uni.json>. In the main directory there will be a statistics hash in JSON format called C<blast_ref.counts.json>.

=head2 Special Methods

=head3 new

    my $blast = Bin::Blast->new($shrub, $workDir, \@refGenomes, \@roles, %options);

Create a new BLAST database in the specified working directory. If the database already exists, it will be reused.

=over 4

=item shrub

The L<Shrub> object for accessing the database.

=item workDir

Name of the working directory to contain the database.

=item refGenomes

Reference to a list of reference genome IDs.

=item roles

Reference to a list of universal role IDs.

=item options

A hash of options. These include

=over 8

=item minsim

The minimum acceptable percent identity similarity. The default is C<0.9>.

=item minlen

The minimum acceptable match length. The default is C<150>.

=item force

If TRUE, the blast database will be rebuilt even if it is already found in the work directory.

=item priv

Privilege level for the functions of the feature. The default is C<1>.

=back

=back

=cut

sub new {
    my ($class, $shrub, $workDir, $refGenomes, $roles, %options) = @_;
    # Extract the options.
    my $minsim = $options{minsim} // 0.9;
    my $minlen = $options{minlen} // 150;
    my $priv = $options{priv} // 1;
    # Verify the work directory.
    if (! $workDir) {
        die "Missing work directory for BLAST database.";
    } elsif (! -d $workDir) {
        die "Invalid BLAST work directory $workDir.";
    }
    # Get the force option.
    my $force = $options{force};
    # Build a hash of universal role IDs.
    my %uniRoles = map { $_ => 1 } @$roles;
    # This will contain the statistics hash.
    my $stats = {};
    # Check to see if we already have statistics. Note that in force mode we start back at zero.
    my $statsFile = "$workDir/blast_ref.counts.json";
    if (-f $statsFile && ! $force) {
        $stats = SeedUtils::read_encoded_object($statsFile);
    }
    # Save the universal role count.
    $stats->{uniRoles} = scalar @$roles;
    # Loop through the reference genomes.
    for my $genome (@$refGenomes) {
        $stats->{refGenomes}++;
        # Get this genome's sub-directory and compute the file names.
        my $genomeDir = "$workDir/$genome";
        my $fastaFile = "$genomeDir/blast_ref.fa";
        my $uniFile = "$genomeDir/blast_ref.uni.json";
        # Check to see if we need to build (or rebuild) the FASTA file.
        if (! -d $genomeDir) {
            # Here the directory doesn't even exist.
            mkdir($genomeDir, 0777);
            build_genome_db($shrub, $genomeDir, $genome, $fastaFile, $uniFile, $statsFile, \%uniRoles, $stats, $priv);
        } elsif ($force || ! -f $fastaFile || ! -f $uniFile) {
            # We have a directory, but we still need a database.
            build_genome_db($shrub, $genomeDir, $genome, $fastaFile, $uniFile, $statsFile, \%uniRoles, $stats, $priv);
        }
    }
    # Create the final object.
    my $retVal = {
        workDir => $workDir,
        minsim => $minsim,
        minlen => $minlen,
        stats => $stats
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 Query Methods

=head3 stats

    my $stats = $blast->stats;

Return a statistics object containing the blast database statistics.

=cut

sub stats {
    my ($self) = @_;
    my $retVal = Stats->new();
    my $stats = $self->{stats};
    for my $stat (keys %$stats) {
        $retVal->Add($stat => $stats->{$stat});
    }
    return $retVal;
}


=head2 Public Methods

=head3 Process

    my ($genomesL, $rolesL) = $blast->Process($seq, \@refs);

Process the specified contig and determine the closest reference genomes and the list of universal roles.

=over 4

=item seq

Sequence of the contig to process.

=item refs

Reference to a list of reference genomes to process. This must be a subset of the reference genomes in this BLAST
database.

=item RETURN

Returns a 2-element list consisting of (0) a reference to a list of the closest genomes (maximum of 2) and (1) a reference to a list
of the IDs of the universal roles found.

=back

=cut

sub Process {
    my ($self, $seq, $refs) = @_;
    # This will track the best match in each genome. It will map each genome ID to a score.
    my %genomes;
    # This will track the universal roles found.
    my %uniRoles;
    # Loop through the given list of genomes.
    for my $genome (@$refs) {
        # Get the genome directory and file names.
        my $genomeDir = "$self->{workDir}/$genome";
        my $fastaFile = "$genomeDir/blast_ref.fa";
        my $uniFile = "$genomeDir/blast_ref.uni.json";
        # Read the universal role hash.
        my $fids2Uni = SeedUtils::read_encoded_object($uniFile);
        # Compute the blast database name.
        my $dbName = gjo::BlastInterface::get_db($fastaFile, 'blastx', $genomeDir);
        # Get the list of matches from the blastx.
        my $matches = gjo::BlastInterface::blast(['contig', '', $seq], $dbName, 'blastx',
                { minIden => $self->{minsim}, minLen => $self->{minlen}, maxE => 1e-50 });
        if (@$matches) {
            print scalar(@$matches) . " hits found in $genome.\n";
            # Loop through the matches. Note they are pre-filtered by the blast method.
            for my $match (@$matches) {
                # Get the ID of the matched protein and the percent identity. Note we don't have to worry about the denominator
                # of the percent identity fraction being 0, since it is pre-filtered on match length.
                my $fid = $match->id2;
                my $score = $match->iden;
                my $role = $fids2Uni->{$fid};
                # If we have a universal role, remember it.
                if ($role) {
                    $uniRoles{$role}++;
                }
                # If this is the best match for this genome, remember it.
                if ($score > ($genomes{$genome} // 0)) {
                    $genomes{$genome} = $score;
                }
            }
        }
    }
    # Now get the best two genomes. If the third or subsequent genome has the same score as the second, we keep it anyway.
    my @sorted = sort { $genomes{$b} <=> $genomes{$a} } keys %genomes;
    my $maxN = scalar @sorted;
    my @genomesL;
    if ($maxN <= 2) {
        @genomesL = @sorted;
    } else {
        @genomesL = shift @sorted;
        my $g2 = shift @sorted;
        push @genomesL, $g2;
        my $g2score = $genomes{$g2};
        while (scalar(@sorted) && $genomes{$sorted[0]} == $g2score) {
            push @genomesL, shift @sorted;
        }
    }
    # Get the list of universal role IDs.
    my @rolesL = keys %uniRoles;
    # Return the information for this contig.
    return (\@genomesL, \@rolesL);
}

=head2 Internal Utility Methods

=head3 build_genome_db

    build_genome_db($shrub, $genomeDir, $genome, $fastaFile, $uniFile, $statsFile, $uniRoles, $stats, $priv);

Create (or re-create) the BLAST database for a single genome.

=over 4

=item shrub

L<Shrub> object for accessing the database.

=item genomeDir

Directory into which the database should be created.

=item genome

ID of the target genome.

=item fastaFile

Name for the output FASTA file of the genome's proteins.

=item uniFile

Name for the output JSON file of the universal role map.

=item statsFile

Name of the statistics file for the whole database.

=item uniRoles

Hash keyed on universal role ID.

=item stats

Hash used to track the statistics.

=item priv

Privilege level for the feature functional assignments.

=back

=cut

sub build_genome_db {
    my ($shrub, $genomeDir, $genome, $fastaFile, $uniFile, $statsFile, $uniRoles, $stats, $priv) = @_;
    print "Processing BLAST data for $genome.\n";
    # This will be the map of feature IDs to universal roles.
    my $fids2Uni = {};
    # We need to writing the proteins to the FASTA file. If a protein has a universal role, we save it in the hash.
    open(my $oh, ">", $fastaFile) || die "Could not open FASTA output file: $!";
    my @prots = $shrub->GetAll('Feature Protein AND Feature Feature2Function', 'Feature(id) LIKE ? AND Feature2Function(security) = ?',
            ["fig|$genome.peg.%", $priv], 'Feature(id) Protein(sequence) Feature2Function(to-link)');
    # Loop through the proteins found.
    for my $prot (@prots) {
        my ($id, $seq, $role) = @$prot;
        $stats->{fids}++;
        # Output the FASTA.
        my @chunks = unpack("(A64)*", $seq);
        print $oh ">$id $role\n";
        print $oh map { "$_\n" } @chunks;
        # Record the role if it is universal.
        if ($uniRoles->{$role}) {
            $fids2Uni->{$id} = $role;
            $stats->{uniFids}++;
        }
    }
    # Close the output file.
    close $oh;
    # Create the BLAST database.
    gjo::BlastInterface::get_db($fastaFile, 'blastx', $genomeDir);
    # Write out the universal role hash.
    SeedUtils::write_encoded_object($fids2Uni, $uniFile);
    # Checkpoint the statistics.
    SeedUtils::write_encoded_object($stats, $statsFile);
}

1;