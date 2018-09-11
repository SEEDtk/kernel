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
use KmerDb;
use Contigs;
use Stats;

=head1 Build DNA KMER Database for Seed Proteins

    seed_kmer_db.pl [ options ] outFile

This script will build a DNA L<KmerDb> for the seed proteins used for binning.  Each protein will be a group unto itself.
The seed protein FASTA file is protein-based, but the sequence IDs refer to features in L<Shrub> which can be used to extract
the DNA.

=head2 Parameters

The positional parameter is the name of the output file. This will contain the L<KmerDb> in json form.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item kmerSize

Number of base pairs per kmer. The default is C<18>.

=item maxFound

The maximum number of occurrences for a kmer to be considered common. The default is C<10>.

=item seedFile

The name of the seed protein file. The default is C<seedprot.fa> in the global data directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outFile',
        Shrub::script_options(),
        ['kmerSize|kmersize|kmer|K=i', 'kmer size in base pairs', { default => 18 }],
        ['maxFound|maxfound|max|m=i',  'maximum number of kmer occurrences', { default => 10 }],
        ['seedFile|seedfile|f=s',      'seed protein FASTA file', { default => "$FIG_Config::global/seedprot.fa" }]
        );
# Check the parameters.
my ($outFile) = @ARGV;
if (! $outFile) {
    die "No output file specified.";
}
my $stats = Stats->new();
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Create the Kmer database.
print "Creating Kmer database.\n";
my $kmerDb = KmerDb->new(kmerSize => $opt->kmersize, maxFound => $opt->maxfound);
# This counter will be used to form group names.
my $groupNo = 0;
# Open the input file.
my $inFile = $opt->seedfile;
print "Reading from $inFile.\n";
open(my $ih, "<$inFile") || die "Could not open $inFile: $!";
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /^>(\S+)/) {
        my $fid = $1;
        $groupNo++;
        $stats->Add(protIn => 1);
        print "Searching for DNA of $fid.\n";
        # This object will contain the DNA for each contig.
        my $contigs = Contigs->new([]);
        # Get the locations for this feature.
        my @locs = $shrub->fid_locs($fid);
        # Loop through the locations, reading contigs.
        for my $loc (@locs) {
            $stats->Add(locIn => 1);
            my $contigID = $loc->Contig;
            if (! $contigs->present($contigID)) {
                my $contigDna = $shrub->contig_seek($contigID);
                $contigs->AddContig($contigID, '', $contigDna);
                $stats->Add(contigIn => 1);
            }
        }
        # Get the DNA.
        my $dna = $contigs->dna(@locs);
        $kmerDb->AddSequence($fid, $dna, "prot$groupNo");
    }
}
# Finalize and save the kmer database.
print "Finalizing kmer database.\n";
$kmerDb->Finalize();
print "Writing to $outFile.\n";
$kmerDb->Save($outFile);
print "All done.\n" . $stats->Show();