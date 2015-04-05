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
use ClosestCoreSEED;

=head1 Get closest genomes from a calibration set

    get_closest_representative_genomes -g GenomesToCalibrate [ options ] < requests

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
                        [ 'calibrationgenomes|g=s', 'file of calibration genomes' ]
    );
my $ih = ScriptUtils::IH( $opt->input );
my $genomes = $opt->calibrationgenomes;

my %genomeH = map { ($_ =~ /^(\d+\.\d+)/) ? ($1 => 1) : () } `cat $genomes`;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my @functions = <DATA>;
chomp @functions;
my $kmer_hash = &ClosestCoreSEED::load_kmers( \@functions, $shrub, $k, \%genomeH );
while ( my $contigs = &read_contigs )
{
    my $response = &ClosestCoreSEED::process_contigs($contigs,$kmer_hash,$k);
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
Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
Phenylalanyl-tRNA synthetase beta chain (EC 6.1.1.20)
Preprotein translocase secY subunit (TC 3.A.5.1.1)
GTP-binding and nucleic acid-binding protein YchF
Translation initiation factor 2
Signal recognition particle, subunit Ffh SRP54 (TC 3.A.5.1.1)
Histidyl-tRNA synthetase (EC 6.1.1.21)
Methionyl-tRNA synthetase (EC 6.1.1.10)
Isoleucyl-tRNA synthetase (EC 6.1.1.5)
Valyl-tRNA synthetase (EC 6.1.1.9)
Seryl-tRNA synthetase (EC 6.1.1.11)
Alanyl-tRNA synthetase (EC 6.1.1.7)
