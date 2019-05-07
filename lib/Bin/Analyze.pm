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


package Bin::Analyze;

    use strict;
    use warnings;
    use Stats;
    use Bin::Score;

=head1 Community Bin Analysis Object.

This object computes a quality score for a set of bins computed by one of the binning algorithms. A bin will be considered
of good quality if it has a specified minimum number of the universal roles and a specified maximum number of duplicate roles.
The quality score for a bin is 1 if it is a good bin and zero if it is not, plus the number of non-duplicate universal roles
divided by the total number of universal roles plus 0.5 if it is a big bin..

This object has the following fields.

=over 4

=item minUnis

The minimum number of universal roles necessary to be considered a good bin.

=item minLen

The minimum number of base pairs required for a bin to be considered big.

=item totUnis

The total number of universal roles. The default is C<101>.

=item genomes

Refererence to a hash mapping genome IDs to genome names.

=back

=head2 Special Methods

=head3 new

    my $analyzer = Bin::Analyze->new(%options);

Construct a new analysis object.

=over 4

=item options

Hash of tuning options.

=over 8

=item minUnis

Minimum number of universal roles necessary to be considered a good bin. The default is C<80>.

=item totUnis

The total number of universal roles.

=item minLen

The minimum number of base pairs required for a bin to be considered big.

=back

=back

=cut

sub new {
    my ($class, %options) = @_;
    # Get the options.
    my $minUnis = $options{minUnis} // 80;
    my $totUnis = $options{totUnis} // 101;
    my $minLen = $options{minLen} // 500000;
    # Create the analysis object.
    my $retVal = {
        minUnis => $minUnis,
        totUnis => $totUnis,
        minLen => $minLen,
        genomes => {}
    };
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head3 Quality

    my $score = Bin::Analyze::Quality(\@bins, %options);

Analyze a list of bins to determine a quality score. This is a non-object-oriented version that can be used for cases
where only one list of bins is being analyzed.

=over 4

=item bins

Reference to a list of L<Bin> objects.

=item options

Hash of tuning options.

=over 8

=item minUnis

Minimum number of universal roles necessary to be considered a good bin. The default is C<51>.

=back

=item RETURN

Returns a value from indicating the number of quality bins.

=back

=cut

sub Quality {
    my ($bins, %options) = @_;
    my $analyze = Bin::Analyze->new(%options);
    my $retVal = $analyze->Analyze($bins);
    return $retVal;
}


=head3 Report

    my $stats = Bin::Analyze::Report(\@bins);

Produce a statistical report about a list of bins. This will show the number of contigs without any BLAST hits,
the number without universal roles, and the distribution of contig lengths, among other things.

=over 4

=item bins

Reference to a list of L<Bin> objects.

=item RETURN

Returns a L<Stats> object containing useful information about the bins.

=back

=cut

sub Report {
    my ($bins) = @_;
    # Get an analysis object.
    my $analyze = Bin::Analyze->new();
    # Create the return object.
    my $stats =$analyze->Stats($bins);
    # Return the statistics.
    return $stats;
}


=head2 Public Methods

=head3 Stats

    my $stats = $analyzer->Stats(\@bins);

Produce a statistical report about a list of bins. This will show the number of contigs without any BLAST hits,
the number without universal roles, and the distribution of contig lengths, among other things.

=over 4

=item bins

Reference to a list of L<Bin> objects.

=item RETURN

Returns a L<Stats> object containing useful information about the bins.

=back

=cut

sub Stats {
    my ($self, $bins) = @_;
    # Create the return object.
    my $stats = Stats->new('goodBin');
    # Loop through the bins.
    for my $bin (@$bins) {
        # Categorize the size, more or less logarithmically. So, we have a category for each
        # multiple of a million for the megabase-order bins, then one for each multiple of 100K,
        # and so forth.
        my $len = $bin->len;
        my $lenCat;
        if ($len < 1000) {
            $lenCat = '0000K';
        } else {
            my $lenThing = 1000000;
            my $zeroes = "";
            my $xes = "XXXK";
            while ($len < $lenThing) {
                $lenThing /= 10;
                $zeroes .= "0";
                $xes = substr($xes, 1);
            }
            my $cat = int($len / $lenThing);
            $lenCat = "$zeroes$cat$xes";
        }
        $stats->Add("binSize-$lenCat" => 1);
        $stats->Add(letters => $len);
        $stats->Add(bins => 1);
        # Check for no proteins and no blast hits.
        my $genomeCount = scalar $bin->refGenomes;
        if (! $genomeCount) {
            $stats->Add(noBlastHits => 1);
            $stats->Add(noUniProts => 1);
        } else {
            $stats->Add(someBlastHits => 1);
            $stats->Add("blastHits-$lenCat" => 1);
            $stats->Add(refHits => $genomeCount);
            my $uniH = $bin->uniProts;
            my $uniCount = scalar keys %$uniH;
            if (! $uniCount) {
                $stats->Add(noUniProts => 1);
            } else {
                $stats->Add(someUniProts => 1);
                my $uniCat = int($uniCount / 10) . "X";
                $stats->Add("uniProtsFound$uniCat" => 1);
                $stats->Add("uniProts-$lenCat" => 1);
                $stats->Add(uniHits => $uniCount);
                for my $uni (keys %$uniH) {
                    $stats->Add("uni-$uni" => $uniH->{$uni});
                }
            }
        }
        # Check for a good bin.
        my $quality = $self->AnalyzeBin($bin);
        if ($quality >= 1) {
            $stats->Add(greatBin => 1);
        }
        if ($quality > 0.5) {
            $stats->Add(goodBin => 1)
        } else {
            $stats->Add(notGoodBin => 1);
        }
    }
    # Return the statistics.
    return $stats;
}


=head3 Analyze

    my $score = $analyzer->Analyze(\@bins);

Analyze a list of bins to determine a quality score.

=over 4

=item bins

Reference to a list of L<Bin> objects.

=item RETURN

Returns a value from indicating the quality of the bins.

=back

=cut

sub Analyze {
    my ($self, $bins) = @_;
    # Analyze the individual bins.
    my $retVal = 0;
    for my $bin (@$bins) {
        $retVal += $self->AnalyzeBin($bin);
    }
    # Return the score.
    return $retVal;
}


=head3 AnalyzeBin

    my $flag = $analyze->AnalyzeBin($bin);

Return the quality score for the bin.

=over 4

=item bin

L<Bin> object to check for sufficient universal roles.

=item RETURN

Returns the number of non-duplicate universal roles divided by the total number of universal roles, plus 1 if the bin is good,
plus 0.5 if the bin is big.

=back

=cut

sub AnalyzeBin {
    my ($self, $bin) = @_;
    # This will be the return value.
    my $retVal = 0;
    # Get this bin's universal role hash.
    my $uniRoles = $bin->uniProts;
    my $uniCount = scalar(keys %$uniRoles);
    # Count the number of duplicates.
    my $dups = 0;
    for my $uniRole (keys %$uniRoles) {
        if ($uniRoles->{$uniRole} > 1) {
            $dups++;
        }
    }
    # Check the universal role count.
    if ($uniCount >= $self->{minUnis}) {
        $retVal = 1;
    }
    # Check the length.
    if ($bin->len >= $self->{minLen}) {
        $retVal += 0.5;
    }
    # Add the full score.
    $retVal += ($uniCount - $dups) / $self->{totUnis};
    # Return the determination indicator.
    return $retVal;
}

=head3 BinReport

    $analyzer->BinReport($oh, $uniRoles, $binList);

Write a detailed report about the bins. Information about the content of the larger bins will be presented, along
with the standard statistical report from L</Report>.

=over 4

=item oh

Open handle for the output file.

=item uniRoles

Reference to a hash mapping each universal role ID to its description. If undefined, the universal roles will
be read from the database.

=item binList

Reference to a list of L<Bin> objects for which a report is desired.

=back

=cut

sub BinReport {
    my ($self, $oh, $uniRoles, $binList) = @_;
    # This will be a hash mapping each universal role to a hash of the bins it appears in. The bins will be
    # identified by an ID number we assign.
    my %uniBins;
    my $binID = 0;
    # Loop through the bins.
    for my $bin (@$binList) {
        # Compute the bin ID.
        $binID++;
        my $bigBin = $self->BinHeader($bin, $oh, $binID);
        if ($bigBin) {
            # Finally, the universal role list. This hash helps us find the missing ones.
            print $oh "    Universal Roles\n";
            print $oh "    ---------------\n";
            my %unisFound = map { $_ => 0 } keys %$uniRoles;
            my $uniFoundCount = 0;
            my $uniMissingCount = 0;
            my $uniDuplCount = 0;
            # Get the universal role hash for the bin.
            my $binUnis = $bin->uniProts;
            for my $uni (sort keys %$binUnis) {
                my $count = $binUnis->{$uni};
                if ($count) {
                    print $oh "    $uni\t$uniRoles->{$uni}\t$count\n";
                    $unisFound{$uni} = 1;
                    $uniBins{$uni}{$binID} = $count;
                    $uniFoundCount++;
                    if ($count > 1) {
                        $uniDuplCount++;
                    }
                }
            }
            # Now the roles not found.
            if (! $uniFoundCount) {
                print $oh "    NONE FOUND\n";
            } else {
                print $oh "    ---------------\n";
                for my $uni (sort keys %unisFound) {
                    if (! $unisFound{$uni}) {
                        print $oh "    $uni\t$uniRoles->{$uni}\tmissing\n";
                        $uniMissingCount++;
                    }
                }
                print $oh "    ---------------\n";
                print $oh "    $uniFoundCount present, $uniMissingCount missing, $uniDuplCount duplicated.\n";
            }
        }
    }
    # Now output the universal role matrix.
    print $oh "\nUNIVERSAL ROLE MATRIX\n";
    print $oh join("\t", 'Role', map { "bin$_" } (1 .. $binID)) . "\n";
    for my $uni (sort keys %uniBins) {
        print $oh join("\t", $uni, map { $uniBins{$uni}{$_} // ' ' } (1 .. $binID)) . "\n";
    }
    print $oh "\n\n";
    # Finally, the bin statistics.
    my $stats = $self->Stats($binList);
    print $oh "FINAL REPORT\n\n" . $stats->Show();
}

=head3 ContigReport

    $analyzer->ContigReport($oh, $id, $bin, \%contigBins);

Write a report on the contig composition of the specified bin to the
specified output stream. For each contig in the bin, the report will list
its universal roles, closest genome, and coverage distance.

=over 4

=item oh

Open output handle to which the report will be written.

=item id

ID number to display for this bin.

=item bin

L<Bin> object for the assembled bin.

=item contigBins

Reference to a hash mapping each contig ID to a single-contig L<Bin> object for it.

=back

=cut

sub ContigReport {
    # Get the parameters.
    my ($self, $oh, $id, $bin, $contigBins) = @_;
    # Display the header.
    my $bigBin = $self->BinHeader($bin, $oh, $id);
    # Only proceed if this is a big bin.
    if ($bigBin) {
        # Loop through the contigs, computing distances.
        my %scoreV;
        my @contigs = $bin->contigs;
        for my $contig (@contigs) {
            my $contigBin = $contigBins->{$contig};
            $scoreV{$contig} = Bin::Score::Vector($bin, $contigBin);
        }
        # Sort the contigs by score.
        my @sorted = sort { Bin::Score::Cmp($scoreV{$a}, $scoreV{$b}) } @contigs;
        print "Contig\tcoverage\tcscore\troles\n";
        for my $contig (@sorted) {
            my $score = $scoreV{$contig};
            my $contigBin = $contigBins->{$contig};
            print $oh join("\t", $contig, $contigBin->meanCoverage, $score->[0], sort keys %{$contigBin->uniProts}) . "\n";
        }
    }
    print "\n\n";
}

=head3 SetGenomes

    $analyzer->SetGenomes(\%genomes);

Specify a hash containing the IDs and names of all the reference genomes.

=over 4

=item genomes

Reference to a hash mapping each genome ID to its name.

=back

=cut

sub SetGenomes {
    my ($self, $genomes) = @_;
    $self->{genomes} = $genomes;
}


=head2 Internal Methods

=head3 BinHeader

    my $bigBin = $analyzer->BinHeader($bin, $oh, $id);

Write the header for a bin in a bin report. Return TRUE if the bin is big
enough to warrant a more detailed analysis.

=over 4

=item bin

L<Bin> object for which a report header is to be produced.

=item oh

Open handle for the output stream.

=item id

ID number to give to the bin.

=item RETURN

Returns TRUE if the bin is big enough for a more detailed report, else FALSE.

=back

=cut

sub BinHeader {
    my ($self, $bin, $oh, $id) = @_;
    # This will be set to TRUE for a big bin.
    my $retVal;
    # Compute the quality score.
    my $quality = $self->AnalyzeBin($bin);
    # Start the header.
    print $oh "\nBIN $id (from " . $bin->contig1 . ", " . $bin->contigCount . " contigs, " . $bin->len . " base pairs, quality $quality)\n";
    if ($bin->name) {
        my $name = $bin->name;
        my $taxonID = $bin->taxonID;
        print $oh "    NCBI taxon $taxonID: $name\n";
    }
    # Only do a detail report if the bin has multiple contigs.
    if ($bin->contigCount > 1) {
        # List the close reference genomes.
        my @genomes = $bin->refGenomes;
        for my $genome (@genomes) {
            print $oh "    $genome: " . $self->gName($genome) . "\n";
        }
        # Compute the average coverage.
        my $avg = $bin->coverage;
        print $oh "*** Mean coverage is $avg.\n";
        # Determine if the bin is big.
        if ($bin->len >= $self->{minLen}) {
            $retVal = 1;
        }
    }
    return $retVal;
}

=head3 gName

    my $genomeName = $analyzer->gName($genome);

Return the name of the specified genome. The name will be cached for subsequent calls.

=over 4

=item genome

ID of the genome whose name is desired.

=item RETURN

Returns the name of the specified genome.

=back

=cut

sub gName {
    my ($self, $genome) = @_;
    # Get the genome hash.
    my $genomeH = $self->{genomes};
    # Look for the name.
    my $retVal = $genomeH->{$genome};
    if (! $retVal) {
        # Create a fake name.
        $retVal = "Unknown species.";
    }
    # Return the name.
    return $retVal;
}

1;