#!/usr/bin/env perl

=head1 Postprocess a set of BLAST hits against a genome database to generate a set of organism counts.

     compile_org_counts < unique_reads.blastn > reads_per_genome
or
     compile_org_counts --input unique_reads.blastn  --output reads_per_genome

=head2 Parameters

There are no positional parameters.


=head2 Options

=over 4

=item input (i)

Name of input file (if not STDIN).


=item output (o)

Name of output file (if not STDOUT).


=item help (h)

Print usage-message and exit.

=back

=cut


#
#  compile_org_counts < unique_reads.blastn > reads_per_genome
#
#  The structure returned by next_blast_query:
#
#     [ qid, qdef, qlen, [ [ sid, sdef, slen, [ hsp_data, hsp_data, ... ] ],
#                          [ sid, sdef, slen, [ hsp_data, hsp_data, ... ] ],
#                          ...
#                        ]
#     ]
#
#     hsp_data:
#
#     [ scr, exp, p_n, pval, nmat, nid, nsim, ngap, dir, q1, q2, qseq, s1, s2, sseq ]
#        0    1    2    3     4     5    6     7     8   9   10   11   12  13   14
#

use strict;
use P3Utils;
use gjoparseblast;
use Data::Dumper;

my $opt = P3Utils::script_opts(
    '',
    [ 'input|i', 'Input file (if not STDIN)', ],
    [ 'output|o', 'Output file (if not STDOUT)', ],
    );

my $input_file  = $opt->input();
my $output_file = $opt->output();

my $input_fh = \*STDIN;
if ($input_file) {
    open($input_fh, q(<), $input_file)
	or die "Could not read-open '$input_file'";
}

my $output_fh = \*STDOUT;
if ($output_file) {
    open($output_fh, q(>), $output_file)
	or die "Could not write-open '$output_file'";
}

my %cnt;
my $reads;
my $subj_data;
my %org;
my $gid;
my $org;

while ( defined( $_ = gjoparseblast::next_blast_query($input_fh) ) )
{
    #  Each query sequence is marked with the number of reads
    ( $reads ) = $_->[1] =~ /^count = (\d+)/;
    $subj_data = $_->[3]->[0];                        #  First hit
    $gid  = $subj_data->[0];                          #  gid is sid
    ( $org = $subj_data->[1] ) =~ s/ SSU rRNA.*$//;   #  Organism is in sdef
    $cnt{ $gid }  += $reads;
    $org{ $gid } ||= $org;
}

my @data = sort { $b->[2] <=> $a->[2] || lc $a->[1] cmp lc $b->[1] }
           map  { [ $_, $org{$_}, $cnt{$_} ] }
           keys %cnt;

foreach ( @data )
{
    printf $output_fh "%8d%13s  %s\n", @$_[2,0,1];
}
