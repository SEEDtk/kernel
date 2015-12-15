#!/usr/bin/env perl
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

use strict;
use warnings;
use Shrub;
use ScriptUtils;
use Data::Dumper;
use gjoseqlib;
use ClosestReps;
use File::Slurp;
use SeedUtils;

=head1 Get closest genomes from a Set of Representatives

    get_closest_representative_genomes -g RepresentativeSet [ options ] < requests

guts of a server that processes requests for closest representative genomes
to new genomes.  The pool of representative (calibration) genomes comes
from a file controlled by the -g parameter.  Each line must begin with a
genome ID.

The input will be requests of the form:

    contigs-in-fasta
    //
    next-set-of-contigs-in-fasta
    //
    .
    .
    .

=head2 Parameters

The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

-g RepresentativeGenomesFile

=cut

my $k = 10;


# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '',
                     Shrub::script_options(), ScriptUtils::ih_options(),
                        [ 'repgenomes|g=s', 'file of representative genomes' ],
                        [ 'closeData|d=s', 'JSON file of kmers for representative genomes']
    );
my $ih = ScriptUtils::IH( $opt->input );
my ($genomeH, $kmer_hash);
if ($opt->closedata) {
    my $closeData = SeedUtils::read_encoded_object($opt->closedata);
    $kmer_hash = $closeData->{kmers};
    $genomeH = $closeData->{genomes};
} elsif ($opt->repgenomes) {
    my $genomes = $opt->repgenomes;
    $genomeH = { map { ($_ =~ /^(\d+\.\d+)\s+(.+)/) ? ($1 => $2) : () } File::Slurp::read_file($genomes) };
    # Connect to the database.
    my $shrub = Shrub->new_for_script($opt);

    my @functions = <DATA>;
    chomp @functions;
    $kmer_hash = &ClosestReps::load_kmers( \@functions,  $genomeH, $k,  $shrub);
} else {
    die "Must specify either a list of calibration genomes or a closeData json file.";
}
while ( my $contigs = &read_contigs )
{
    my $response = &ClosestReps::process_contigs($contigs,$kmer_hash,$k, $genomeH);
    print $response, "\n//\n";
}

sub read_contigs
{
    my $contigs = [];
    $/ = "\n//\n";
    my $req = <$ih>;
    if ($req)
    {
        chomp $req;
        my @entries = split( /\n\>\s*/, $req );
        foreach my $entry (@entries)
        {
            if ( $entry =~ /^(>\s*)?(\S+)[^\n]*\n(.*)/s )
            {
                my ( $id, $seq ) = ( $2, $3 );
                $seq =~ s/\s+//gs;
                push( @$contigs, [ $id, '', $seq ] );
            }
        }
    }
    $/ = "\n";
    if (@$contigs)
    {
        return $contigs;
    }
    return undef;
}
__DATA__
GtpBindNuclAcidBind
PhenTrnaSyntAlph
PhenTrnaSyntBeta
PrepTranSecySubu
TranInitFact2n1
HistTrnaSynt
SignRecoPartSubu
MethTrnaSynt
IsolTrnaSynt
SeryTrnaSynt
ValyTrnaSynt
AlanTrnaSynt
