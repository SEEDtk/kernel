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
use FIG_Config;
use Shrub;
use ScriptUtils;
use Data::Dumper;
use FIG_Config;
use ClosestCoreSEED;
use gjo::seqlib;

my $k = 10;

=head1 Test the ability of get_closest_coreSEED_genomes.pl

    get_closest_genomes_all.pl > best.picks.from.coreSEED

=head2 Parameters


The command-line options are those found in L<Shrub/script_options> and
L<ScriptUtils/ih_options> plus the following.

=cut

# Get the command-line parameters.
my $opt =
  ScriptUtils::Opts( '', Shrub::script_options(), ScriptUtils::ih_options(),
    [] );
my $ih = ScriptUtils::IH( $opt->input );

# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my $min_hits = 1000;  ### minimum number of kmer hits

my @tuples = $shrub->GetAll("Genome","",[],"Genome(id) Genome(name) Genome(contig-file)");
my $dnaRepo = $shrub->DNArepo();

my @functions = <DATA>;
chomp @functions;
my $kmer_hash = &ClosestCoreSEED::load_kmers( \@functions, $shrub, $k );

foreach my $tuple (@tuples)
{
    my($g,$gs,$cf) = @$tuple;
    my $contig_file = "$dnaRepo/$cf";
    my @contigs = &read_fasta($contig_file);
    my $response = &ClosestCoreSEED::process_contigs(\@contigs,$kmer_hash,$k);
    my @out = split(/\n/,$response);
    my @hits = map { (($_ =~ /^(\d+), (\S+)/) && ($1 >= $min_hits)) ? [$1,$2] : () } @out;
    if (@hits == 0)
    {
	print STDERR join("\t",('could not be placed',$g,$gs)),"\n";
    }
    else
    {
	print join("\t",($hits[0]->[0],$hits[0]->[1],$g,$gs)),"\n";
	foreach my $hit (@hits)
	{
	    print join("\t",@$hit),"\n";
	}
	print "//\n";
    }
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
