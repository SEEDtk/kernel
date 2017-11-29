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


package BinningReports;

    use strict;
    use warnings;
    use URI::Escape;
    use File::Basename;
    use Template;
    use Data::Dumper;

=head1 Produce Binning Reports

This package produces binning reports from the output JSON objects. The main method takes as input the various objects and the two templates and
produces an HTML string for the summary report plus one for each detail report. The method can be tested in a batch testing environment and then
plugged into the main PATRIC binning engine.

=head2 Public Methods

=head3 Process

    my ($summary, $detailsH) = BinningReports::Process($jobID, $params, $quality_json, $ppr_json, $bins_json, $summary_tt, $detail_tt, $genome_group_path);

Produce the required quality reports. This includes one summary report plus a detail report for each bin.

=over 4

=item jobID

The PATRIC job identifier for the binning run.

=item params

The parameter structure used to invoke the binning. This should be a hash reference with the following keys. All are optional; however, at least
one of C<contigs> and C<paired_end_libs> must be present.

=over 8

=item contigs

Name of the contigs input file.

=item paired_end_libs

Reference to a list containing the two paired-end read input files.

=item genome_group

Name of the genome group into which output genomes are to be placed.

=back

=item quality_json

A reference to a hash containing the following keys.

=over 8

=item package_directory

The path to the output genome packages.

=item packages

A list of hash references, one per bin, containing the following keys.

=over 8

=item checkm_completeness

The CheckM completeness score of the bin's genome, in percent. High values are good.

=item checkm_contamination

The CheckM contamination score of the bin's genome, in percent. High values are bad.

=item contigs

The number of contigs in the bin's genome.

=item genome_id

The ID of the bin's genome in PATRIC. Each genome ID is completely unique.

=item genome_name

The name given to the bin's genome by the binning software. Each bin will have a name unique to the binning run.

=item n50

The N50 measure of contig size distribution of the bin's genome. A higher number is an indication of better assembly.

=item scikit_coarse

The coarse SciKit consistency score of the bin's genome, in percent. This roughly corresponds to completeness.

=item scikit_fine

The fine SciKit consistency score of the bin's genome, in percent. This is always lower than the coarse score, and is influenced
by contamination.

=back

=back

=item ppr_json

A reference to a hash containing that maps each genome ID to a sub-hash containing the following keys.

=over 8

=item roles

Reference to a hash that maps each role ID to its name.

=item role_fids

Reference to a hash that maps each role ID to a list of feature IDs for features possessing that role.

=item role_ppr

Reference to a hash that maps each role ID to a 2-tuple consisting of (0) the number of expected occurrences and (1) the number of actual occurrences.

=back

=item bins_json

Reference to a list of hash references, one per bin, with the following keys.

=over 4

=item binFastaFile

The name of the file containing the contigs for the bin.

=item binFastaPath

The directory path of the file containing the contigs for the bin.

=item binIndex

The index number of the bin (1-based).

=item contigs

Reference to a list of the IDs of the contigs in the bin.

=item coverage

The mean coverage for the bin.

=item domain

The domain of the bin's genome.

=item geneticCode

The majority genetic code used for the bin.

=item len

The number of DNA base pairs in the bin.

=item name

The name assigned to the bin's genome.

=item refGenomes

A reference to a list of the IDs for the bin's reference genomes.

=item taxonID

The estimated taxonomic ID for the bin.

=item tetra

A reference to a vector of tetranucleotide counts for the bin's DNA.

=item tetraLen

The length of the tetranucleotide vector for the bin, used to normalize it.

=item uniProts

Reference to an empty hash.

=back

=item summary_tt

The template to be used for the summary report. This template expects the following variables.

=over 8

=item job_id

The PATRIC job identifier.

=item min_checkm

The minimum CheckM completeness score required for a good bin.

=item min_scikit

The minimum SciKit completeness score required for a good bin.

=item max_contam

The maximum CheckM contamination score allowed for a good bin.

=item params

The I<params> structure described above.

=item found

Reference to a hash with the following keys.

=over 12

=item total

Total number of bins.

=item good

Number of good bins.

=item bad

Number of bad bins.

=back

=item good

Reference to a list of bin descriptors for good bins. Each descriptor contains the bin's entry from the I<quality_json> structure
described above, plus the following members.

=over 12

=item ppr

A count of the problematic roles.

=item refs

A list of reference genome descriptors, each a hash reference containing the genome name in C<genome> and the URL in C<url>.

=item coverage

The mean coverage for the bin.

=back

=item bad

Reference to a list of bin descriptors for bad bins. The descriptors are identical in format to the ones for the good bins.

=back

=item details_tt

The template to be used for each bin's detail report. This template expects the following variables.

=over 8

=item g

A bin descriptor for the bin, as described in the list of good-bin descriptors above.

=item p

A problematic-role descriptor for the bin, consisting of a list of structures, one per role. Each structure contains the following fields.

=over 12

=item role

The role description.

=item predicted

The number of roles predicted.

=item actual

The actual number of roles found.

=item n_fids

The number of features containing the problematic role.

=item fid_url

The URL to list the features.

=back

=back

=item genome_group_path

The path to the output genome group (if any).

=item roleMap

A reference to a hash mapping role IDs to role descriptions.

=item RETURN

Returns a two-element list. The first is the HTML for the summary report. The second is a reference to a hash mapping each bin's genome ID
to the HTML string for its detail report.

=back

=cut

# Good/Bad criteria
use constant MIN_CHECKM => 80;
use constant MIN_SCIKIT => 85;
use constant MAX_CONTAM => 10;

# URL helpers
use constant URL_BASE => 'https://www.patricbrc.org/view/Genome';
use constant FID_URL_BASE => 'https://www.patricbrc.org/view/Genome';

sub Process {
    my ($jobID, $params, $quality_json, $ppr_json, $bins_json, $summary_tt, $detail_tt, $genome_group_path, $roleMap) = @_;
    # We need to create $found, $good, and $bad for the summary report. We also need to create an object for each
    # detail report. The detail report objects go into this hash.
    my %detailParms;
    # Here are the storage places for found, good, and bad. The "g" item from the bin's descriptor goes in these
    # lists.
    my %found = (total => 0, good => 0, bad => 0);
    my (@good, @bad);
    # First we are going to read through the bins and create a map of bin names to reference genome descriptors and coverages.
    # Each reference genome descriptor is a hash-ref with members "genome" and "url".
    my %refGmap;
    for my $binThing (@$bins_json) {
        my $name = $binThing->{name};
        my $refs = $binThing->{refGenomes};
        my @refList = map { { genome => $_, url => join('/', URL_BASE , uri_escape($_)) } } @$refs;
        my ($cov, $count) = (0, 0);
        for my $covItem (@{$binThing->{coverage}}) {
            $cov += $covItem;
            $count++;
        }
        if ($count > 1) {
            $cov /= $count;
        }
        $cov = int($cov * 100) / 100;
        $refGmap{$name} = { refs => \@refList, coverage => $cov };
    }
    # The quality_json is the master bin list. As we process it, we connect each entry to its matching data in ppr_report.json and
    # %refGmap.
    my $packagesL = $quality_json->{packages};
    my @packages = sort { quality_score($b) <=> quality_score($a) } @$packagesL;
    for my $bin (@packages) {
        # Copy the package entry. This copy will be made into the main object used to describe bins in the output reports.
        my %gThing = %$bin;
        # Get the matching ppr and refGmap entries.
        my $genomeID = $bin->{genome_id};
        my $genomeName = $bin->{genome_name};
        my $genomeURL = join('/', URL_BASE, uri_escape($genomeID));
        my $ppr = $ppr_json->{$genomeID};
        my $refData = $refGmap{$genomeName};
        # Only proceed if we connected the pieces. We need a fall-back in case of errors.
        if ($ppr && $refData) {
            # Connect the coverage and reference genome data.
            $gThing{refs} = $refData->{refs};
            $gThing{coverage} = $refData->{coverage};
            $gThing{genome_url} = $genomeURL;
            # Create the ppr descriptor for the bin. This also nets us the ppr count.
            my $pprs = 0;
            my $pprRoleFids = $ppr->{role_fids};
            my $pprRoleData = $ppr->{role_ppr};
            my @pprList;
            for my $role (sort keys %$pprRoleData) {
                my $pa = $pprRoleData->{$role} // [0,0];
                my ($predicted, $actual) = @$pa;
                if ($predicted != $actual) {
                    $pprs++;
                    my $roleName = $roleMap->{$role} // $role;
                    my $fidList = $pprRoleFids->{$role} // [];
                    my $n_fids = scalar @$fidList;
                    my %pprThing = (role => $roleName, predicted => $predicted, actual => $actual, n_fids => $n_fids);
                    $pprThing{fid_url} = fid_list_url($fidList);
                    push @pprList, \%pprThing;
                }
            }
            # Store the PPR count in the main descriptor.
            $gThing{ppr} = $pprs;
            # Attach this data to the bin.
            $detailParms{$genomeID} = { g => \%gThing, p => \@pprList };
            # Is this bin good or bad?
            if ($gThing{checkm_completeness} >= MIN_CHECKM && $gThing{scikit_fine} >= MIN_SCIKIT &&
                    $gThing{checkm_contamination} <= MAX_CONTAM) {
                push @good, \%gThing;
                $found{good}++;
            } else {
                push @bad, \%gThing;
                $found{bad}++;
            }
            # Update the total-bin count.
            $found{total}++;
        }
    }
    # We have now compiled the information we need for each report. Create the template engine.
    my $templateEngine = Template->new(ABSOLUTE => 1);
    # Allocate the result variables.
    my ($summary, %detailsH);
    # Create the summary report parm structure.
    my $vars = { job_id => $jobID, params => $params, found => \%found, good => \@good, bad => \@bad, group_path => $genome_group_path,
                 min_checkm => MIN_CHECKM, min_scikit => MIN_SCIKIT, max_contam => MAX_CONTAM };
    # print STDERR Dumper($vars);
    # Create the summary report.
    $templateEngine->process($summary_tt, $vars, \$summary);
    # Loop through the bins.
    for my $bin (keys %detailParms) {
        $vars = $detailParms{$bin};
        my $detail;
        # print STDERR Dumper($vars);
        $templateEngine->process($detail_tt, $vars, \$detail);
        $detailsH{$bin} = $detail;
    }
    # Return the reports.
    return ($summary, \%detailsH);
}

=head3 quality_score

    my $sortVal = BinningReports::quality_score($g);

Determine the quality score for a package object from the quality_json structure. A higher quality score means a better bin.

=over 4

=item g

The package object for the bin in question.

=item RETURN

Returns a score that is higher for better bins.

=back

=cut

sub quality_score {
    my ($g) = @_;
    my $retVal = $g->{checkm_completeness} + 1.1 * $g->{scikit_fine} - 5 * $g->{checkm_contamination};
    return $retVal;
}

=head3 fid_list_url

    my $url = BinningReports::fid_list_url(\@fids);

Return a URL for viewing a list of PATRIC features.

=over 4

=item fids

Reference to a list of feature IDs.

=item RETURN

Returns a URL for viewing a single feature or a list of multiple features.

=back

=cut

sub fid_list_url {
    my ($fids) = @_;
    my $retVal;
    if (@$fids == 1) {
        $retVal = "https://www.patricbrc.org/view/Feature/" . uri_escape($fids->[0]);
    } elsif (@$fids > 1) {
        my $list = join(",", map { uri_escape(qq("$_")) } @$fids);
        $retVal = "https://www.patricbrc.org/view/FeatureList/?in(patric_id,($list))";
    }
    return $retVal;

}

1;