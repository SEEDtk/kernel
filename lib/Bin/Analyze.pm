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

=head1 Community Bin Analysis Object.

This object computes a quality score for a set of bins computed by the L<Bin::Compute> algorithm. A bin will be considered
of good quality if it has a specified minimum number of the universal roles and a specified maximum number of duplicate roles.
The quality score is the number of good bins.

This object has the following fields.

=over 4

=item minUnis

The minimum number of universal roles necessary to be considered a good bin.

=item maxDups

The maximum number of duplicate universal roles allowed in a good bin.

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

Minimum number of universal roles necessary to be considered a good bin. The default is C<51>.

=item maxDups

Maximum number of duplicate universal roles allowed in a good bin. The default is C<4>.

=back

=back

=cut

sub new {
    my ($class, %options) = @_;
    # Get the options.
    my $minUnis = $options{minUnis} // 30;
    my $maxDups = $options{maxDups} // 4;
    # Create the analysis object.
    my $retVal = {
        minUnis => $minUnis,
        maxDups => $maxDups
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

=item maxDups

Maximum number of duplicate universal roles allowed in a good bin. The default is C<4>.

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
        if ($self->AnalyzeBin($bin)) {
            $stats->Add(goodBin => 1);
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

Returns a value from 0 to 1 indicating the proportion of quality bins.

=back

=cut

sub Analyze {
    my ($self, $bins) = @_;
    # Analyze the individual bins.
    my @good = grep { $self->AnalyzeBin($_) == 2 } @$bins;
    # Return the number of good ones.
    return scalar @good;
}


=head3 AnalyzeBin

    my $flag = $analyze->AnalyzeBin($bin);

Return 2 if the specified bin is good, 1 if it has a lot of universal roles, else 0.

=over 4

=item bin

L<Bin> object to check for sufficient universal roles.

=item RETURN

Returns 2 if the bin has sufficient universal roles and not too many duplicates, 1 if it has sufficient universal roles
and too many duplicates, else 0.

=back

=cut

sub AnalyzeBin {
    my ($self, $bin) = @_;
    # This will be the return value.
    my $retVal = 0;
    # Get this bin's universal role hash.
    my $uniRoles = $bin->uniProts;
    # Check the universal role count.
    if (scalar(keys %$uniRoles) >= $self->{minUnis}) {
        # It's good. Count the number of duplicates.
        my $dups = 0;
        for my $uniRole (keys %$uniRoles) {
            if ($uniRoles->{$uniRole} > 1) {
                $dups++;
            }
        }
        if ($dups <= $self->{maxDups}) {
            $retVal = 2;
        } else {
            $retVal = 1;
        }
    }
    # Return the determination indicator.
    return $retVal;
}


1;