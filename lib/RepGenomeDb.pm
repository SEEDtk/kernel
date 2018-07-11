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


package RepGenomeDb;

    use strict;
    use warnings;
    use RepGenome;
    use gjoseqlib;

=head1 Representative Genome Database

This object encapsulates a set of representative genomes. It includes the basic parameters of the representation (kmer size, minimum similarity score),
one L<RepGenome> object per representative genome, and a hash mapping each represented genome to its representative.

The fields of this object are as follows.

=over 4

=item gMap

Reference to a hash mapping the ID of each reference genome to its L<RepGenome> object.

=item repMap

Reference to a hash mapping each represented genome to its reference genome.

=item K

Kmer size to use, in amino acids.

=item score

The minimum similarity score for two genomes to be considered close.

=back

=head2 Special Methods

=head3 new

    my $repDB = RepGenomeDb->new(%options);

Create a new, blank representative-genome database.

=over 4

=item options

A hash of configuration options, containing zero or more of the following keys.

=over 8

=item K

The kmer size to use for similarity computation. The default is C<8>.

=item score

The minimum score (kmers in common) for two genomes to be considered close. The default is C<100>.

=back

=back

=cut

sub new {
    my ($class, %options) = @_;
    # Extract the options.
    my $k = $options{K} // 8;
    my $score = $options{score} // 100;
    # Form the empty database.
    my $retVal = { K => $k, score => $score, repMap => {}, gMap => {} };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}

=head3 new_from_dir

    my $repDB = RepGenome->new_from_dir($dirName, %options);

Read a representative-genome database from a directory.

=over 4

=item dirName

The name of a representative-genome database directory. The directory must contain the following files.

=over 8

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item rep_db.tbl

A tab-delimited file containing one record per genome with the columns (0) the genome ID, (1) the ID of the genome's
representative, and (2) the similarity number. If this file is missing, it will be presumed empty.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

=item options

A hash containing zero or more of the following keys.

=over 8

=item verbose

If TRUE, status messages will be sent to the standard output.

=back

=back

=cut

sub new_from_dir {
    my ($class, $dirName, %options) = @_;
    # Check the options.
    my $verbose = $options{verbose};
    # Load the parameter file.
    my ($k, $score);
    if (! -s "$dirName/K") {
        # No parm file, so default the parameters.
        ($k, $score) = (8, 100);
    } else {
        # Here we have a parm file.
        open(my $kh, '<', "$dirName/K") || die "Could not open parameter file for $dirName: $!";
        ($k, $score) = map { $_ =~ /(\d+)/; $1 } <$kh>;
        # Older parm files do not have the score, so default it.
        $score //= 100;
    }
    # Create the object.
    my $retVal = new($class, K => $k, score => $score);
    # We need to create the representative-genome map. First we get the ID and name of each one.
    print "Reading genome names.\n" if $verbose;
    open(my $kh, '<', "$dirName/complete.genomes") || die "Could not open genome file for $dirName: $!";
    my $gNames;
    while (! eof $kh) {
        my $line = <$kh>;
        if ($line =~ /^(\d+\.\d+)\s+(.+)/) {
            $gNames->{$1} = $2;
        }
    }
    close $kh; undef $kh;
    # Now we connect the proteins.
    $retVal->_ReadProteins("$dirName/6.1.1.20.fasta", $gNames, $verbose);
    # Finally, we read the rep_db.tbl file to get the represented genomes. If the file doesn't exist we skip this step.
    my $repDbFile = "$dirName/rep_db.tbl";
    if (-s $repDbFile) {
        print "Reading represented-genomes list.\n" if $verbose;
        open($kh, '<', $repDbFile) || die "Could not open represented-genomes file for $dirName: $!";
        while (! eof $kh) {
            my $line = <$kh>;
            my ($genome, $repID, $score) = $line =~ /^(\d+\.\d+)\t(\d+\.\d+)\t(\d+)/;
            if ($score) {
                $retVal->Connect($repID, $genome, $score);
            }
        }
    }
    # Return the object created.
    return $retVal;
}

=head2 Query Methods

=head3 K

    my $K = $self->K();

Return the kmer size for this database.

=cut

sub K {
    my ($self) = @_;
    return $self->{K};
}

=head3 score

    my $score = $self->score();

Return the minimum similarity score for this database. To be represented, a genome must have this score or better with a representative.

=cut

sub score {
    my ($self) = @_;
    return $self->{score};
}

=head3 check_rep

    my ($repID, $score) = $repDB->check_rep($genome);

If the specified genome is currently represented, return the ID of its representative and the similarity score.

=over 4

=item genome

ID of the genome to check.

=item RETURN

Returns a two-element list consisting of (0) the ID of the representative genome, and (1) the similarity score. If the
genome is not currently represented, returns C<undef> for the genome ID and C<0> for the score.

=back

=cut

sub check_rep {
    my ($self, $genome) = @_;
    # Create the result variables.
    my ($repID, $score) = (undef, 0);
    # Get the representing genome's object.
    my $repGenome = $self->{repMap}{$genome};
    if ($repGenome) {
        # We found one, so extract the ID and similarity score.
        $repID = $repGenome->id();
        $score = $repGenome->score($genome);
    }
    # Return the results.
    return ($repID, $score);
}


=head3 find_rep

    my ($repID, $score) = $repDB->find_rep($prot);

Find the representative genome closest to the specified protein and return the similarity score.

=over 4

=item prot

The identifying protein sequence to use for the similarity determination.

=item RETURN

Returns a two-element list consisting of (0) the ID of the closest genome and (1) the similarity score.

=back

=cut

sub find_rep {
    my ($self, $prot) = @_;
    # Initialize the best score and the closest-genome ID.
    my ($repID, $score) = (undef, 0);
    # Loop through the genomes, remembering the best match.
    my $gMap = $self->{gMap};
    for my $genome (keys %$gMap) {
        my $repGenome = $gMap->{$genome};
        my $gScore = $repGenome->check_genome($prot);
        if ($gScore > $score) {
            $repID = $genome;
            $score = $gScore;
        }
    }
    # Return the best match.
    return ($repID, $score);
}

=head3 count_rep

    my $count = $repDB->count_rep($prot, $score);

Count the number of representative genomes for the specified protein with the specified similarity or better.

=over 4

=item prot

The identifying protein sequence to use for the similarity determination.

=item score

The minimum acceptable similarity score.

=item RETURN

Returns the number of representative genomes closer to the protein than the specified score.

=back

=cut

sub count_rep {
    my ($self, $prot, $score) = @_;
    # Initialize the count.
    my $retVal = 0;
    # Loop through the genomes, remembering the best match.
    my $gMap = $self->{gMap};
    for my $genome (keys %$gMap) {
        my $repGenome = $gMap->{$genome};
        my $gScore = $repGenome->check_genome($prot);
        if ($gScore >= $score) {
            $retVal++;
        }
    }
    # Return the count.
    return $retVal;
}


=head3 list_reps

    my $repList = $repDB->list_reps($prot, $score);

List the representative genomes for the specified protein with the specified similarity or better.

=over 4

=item prot

The identifying protein sequence to use for the similarity determination.

=item score

The minimum acceptable similarity score.

=item RETURN

Returns a reference to a hash mapping the IDs of close genomes to their similarity score.

=back

=cut

sub list_reps {
    my ($self, $prot, $score) = @_;
    # This will be the return hash.
    my %retVal;
    # Loop through the genomes, looking for good matches.
    my $gMap = $self->{gMap};
    for my $genome (keys %$gMap) {
        my $repGenome = $gMap->{$genome};
        my $gScore = $repGenome->check_genome($prot);
        if ($gScore >= $score) {
            $retVal{$genome} = $gScore;
        }
    }
    # Return the hash.
    return \%retVal;
}


=head3 rep_list

    my $repHash = $repDB->rep_list();

Return a list of the IDs of the representative genomes in this database.

=cut

sub rep_list {
    my ($self) = @_;
    my $gMap = $self->{gMap};
    return [sort keys %$gMap];
}

=head3 rep_object

    my $repGenome = $repDB->rep_object($genomeID);

Return the L<RepGenome> object for the genome with the specified ID, or C<undef> if the genome is not found.

=over 4

=item genomeID

The ID of the genome whose L<RepGenome> object is desired.

=item RETURN

Returns the L<RepGenome> object for the identified genome, or C<undef> if the genome is not in the database.

=back

=cut

sub rep_object {
    my ($self, $genomeID) = @_;
    my $gMap = $self->{gMap};
    my $retVal = $gMap->{$genomeID};
    return $retVal;
}

=head2 Public Manipulation Methods

=head3 AddRep

    $repDB->AddRep($id, $name, $prot);

Add a representative genome to the database.

=over 4

=item id

The ID of the genome to add.

=item name

The name of the genome.

=item prot

The sequence of the genome's identifying protein.

=back

=cut

sub AddRep {
    my ($self, $id, $name, $prot) = @_;
    # Get the kmer size.
    my $K = $self->{K};
    # Create the representative-genome object.
    my $repGenome = RepGenome->new($id, name => $name, prot => $prot, K => $K);
    # Add it to the rep-map.
    $self->{gMap}{$id} = $repGenome;
}

=head3 Connect

    $repDB->Connect($repID, $genome, $score);

Connect a represented genome to its representative.

=over 4

=item repID

The ID of the representative genome.

=item genome

The ID of the genome being represented.

=item score

The similarity score between the two genomes.

=back

=cut

sub Connect {
    my ($self, $repID, $genome, $score) = @_;
    # Get the two maps.
    my $repMap = $self->{repMap};
    my $gMap = $self->{gMap};
    # Get the representative-genome object for the representative genome.
    my $repGenome = $gMap->{$repID} // die "$repID not found in representative-genome database.";
    # Add this genome to it.
    $repGenome->AddGenome($genome, $score);
    # Denote it is the representative of this genome.
    $repMap->{$genome} = $repGenome;
}

=head3 Save

    $repDB->Save($outDir);

Save a copy of the new database to the specified output directory.

=over 4

=item outDir

The directory into which the input files for re-creating this database will be written.

=back

=cut

sub Save {
    my ($self, $outDir) = @_;
    # First, write the parameter file.
    my $K = $self->{K};
    my $score = $self->{score};
    open(my $oh, '>', "$outDir/K") || die "Could not open output parameter file for $outDir: $!";
    print $oh "$K\n$score\n";
    close $oh; undef $oh;
    # Now we need to run through the gMap. We create three files in parallel-- 6.1.1.20.fasta, rep_db.tbl, and complete.genomes.
    open($oh, '>', "$outDir/complete.genomes") || die "Could not open output genome file for $outDir: $!";
    open(my $fh, '>', "$outDir/6.1.1.20.fasta") || die "Could not open output FASTA file for $outDir: $!";
    open(my $rh, '>', "$outDir/rep_db.tbl") || die "Could not open output rep_db file for $outDir: $!";
    my $gMap = $self->{gMap};
    for my $genome (sort keys %$gMap) {
        my $repGenome = $gMap->{$genome};
        # Put the name in complete.genomes.
        my $name = $repGenome->name();
        print $oh "$genome\t$name\n";
        # Put the identifying protein's sequence in 6.1.1.20.fasta.
        my $prot = $repGenome->prot();
        print $fh ">$genome\n$prot\n";
        # Put the represented genomes in rep_db.tbl.
        my $repList = $repGenome->rep_list();
        for my $repTuple (@$repList) {
            print $rh join("\t", $repTuple->[0], $genome, $repTuple->[1]) . "\n";
        }
    }
    # Close all the files.
    close $oh;
    close $fh;
    close $rh;
}

=head3 AddRepObject

    $repDB->AddRepObject($repGenome);

Add a new representative genome via a L<RepGenome> object.

=over 4

=item repGenome

A L<RepGenome> object for the representative genome to add.

=back

=cut

sub AddRepObject {
    my ($self, $repGenome) = @_;
    # Get the two maps.
    my $repMap = $self->{repMap};
    my $gMap = $self->{gMap};
    # Get the genome ID.
    my $id = $repGenome->id();
    # Add the object to the rep-map.
    $gMap->{$id} = $repGenome;
    # Get its connected genomes.
    my $repList = $repGenome->rep_list();
    # Put them in the cross-reference map.
    for my $tuple (@$repList) {
        my ($genomeID) = $tuple->[0];
        $repMap->{$genomeID} = $repGenome;
    }
}

=head2 Internal Utilities

=head3 _ReadProteins

    $repDB->_ReadProteins($fileName, \%gNames);

Create the representative-genome map from a hash of genome names and the protein FASTA file.

=over 4

=item fileName

Name of the protein FASTA file. Each record should contain a genome ID and the sequence of the identifying protein. The genome
ID is normally taken from the comment field, but in the absence of a comment the sequence ID will be used.

=item gNames

Reference to a hash mapping genome IDs to genome names. Only genomes identified in this hash will be added to the map.

=item verbose

If TRUE, trace messages will be displayed.

=back

=cut

sub _ReadProteins {
    my ($self, $fileName, $gNames, $verbose) = @_;
    # Read the FASTA file.
    open(my $kh, '<', $fileName) || die "Could not open protein FASTA $fileName: $!";
    my @triples = gjoseqlib::read_fasta($kh);
    print scalar(@triples) . " proteins found in FASTA.\n" if $verbose;
    for my $triple (@triples) {
        my ($id, $genome, $seq) = @$triple;
        $genome ||= $id;
        # Now we have a genome ID and a sequence. Get the genome name as well.
        my $name = $gNames->{$genome};
        # If the genome is in the FASTA but we don't have a name, then it's an extra genome. (We allow this as
        # a convenience.) If we have the name, we process it.
        if ($name) {
            $self->AddRep($genome, $name, $seq);
        }
    }
    undef @triples;
}

1;