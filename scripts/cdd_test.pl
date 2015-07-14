use strict;
use Shrub;
use ScriptUtils;
use Data::Dumper;
use Carp;
use Getopt::Long;
use LWP::UserAgent;
use Digest::MD5;
use ConservedDomainSearch;

use gjoseqlib;

my $usage = <<"End_of_Usage";

usage: svr_cdd_scan [options] < seqs.fa > cdd.table

       -h   print this help screen
       -s   print short output table (one line for each sequence):

  Short output:

  [ seq_id, specific_hits, superfamiy_hits, multi-domain_hits, 
            names_for_specific_hits, names_for_superfamiy_hits, names_for_multi-domain_hits ]

  Long output (default):

  [ seq_id query_no query_id hit_type pssm_id from to e_value bit_score
           accession domain short_name incomplete_flag superfamily ] 

End_of_Usage

my ($help, $direct, $short, $idlist);

ScriptUtils::Opts( '',
                      Shrub::script_options(), ScriptUtils::ih_options(),
                                           [ 'short|s', 'produce short output']
                                               );

$help and die $usage;

my %options =  (data_mode => 'rep');

my $md5;

my $CDD = ConservedDomainSearch::new();
print Dumper $CDD; die();

while (my ($id, $def, $seq) = read_next_fasta_seq()) {
        my $d = $CDD->lookup_seq($id, $md5, $seq, \%options);
        print $id, "\n";
        print Dumper $d;
}
