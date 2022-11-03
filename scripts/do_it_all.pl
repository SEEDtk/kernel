#!/usr/bin/perl
#
#  do_it_all  forward_reads  reverse_reads  SSU_rRNA  [basename]
#
use strict;
use Time::HiRes qw( time );

my $example = <<'End_of_example';

do_it_all FIG-Amplicon_8_12_21_S1_L001_R1_001.fna FIG-Amplicon_8_12_21_S1_L001_R2_001.fna rep100_SSU_rRNA_515F-926R tmp2

End_of_example

my $usage = <<'End_of_usage';
End_of_usage

my $keep_blast = 1;

my ( $fwd_read, $rev_read, $SSU_rRNA, $basename ) = splice( @ARGV );

-s $fwd_read
    or print STDERR "Forward reads not found.\n", $usage
        and exit;

-s $rev_read
    or print STDERR "Reverse reads not found.\n", $usage
        and exit;

-s "$SSU_rRNA.nsq"
    or print STDERR "SSU rRNA blast database not found.\n", $usage
        and exit;

if ( ! $basename )
{
    my $i = 1;
    while ( 1 )
    {
        $basename = sprintf "rRNA_analysis_%09d", $i;
        last if ! ( -f "$basename.joined"
                 || -f "$basename.unique"
                 || -f "$basename.blastn"
                 || -f "$basename.org_counts"
                  );
        $i++;
    }
}

my $t0 = time();

print STDERR "\n";
print STDERR "Running overlap_reads ...\n";
system( qq(overlap_reads -s '$fwd_read' '$rev_read' > '$basename.joined') );  # -s stringent overlap
my $t1 = time();
printf STDERR "... %.3f seconds\n\n", $t1-$t0;

print STDERR "Running make_nonredundant ...\n";
system( qq(make_nonredundant < '$basename.joined' > '$basename.unique') );
my $t2 = time();
printf STDERR "... %.3f seconds\n\n", $t2-$t1;

print STDERR "Running blastn ...\n";
system( qq(blastn -query '$basename.unique' -db '$SSU_rRNA' -evalue 1e-40 -perc_identity 40 -qcov_hsp_perc 80 -penalty -1 -reward 1 -gapopen 2 -gapextend 1 -num_descriptions 5 -num_alignments 5 -num_threads 8 > '$basename.blastn') );
my $t3 = time();
printf STDERR "... %.3f seconds\n\n", $t3-$t2;

print STDERR "Running compile_org_counts ...\n";
system( qq(compile_org_counts < '$basename.blastn' > '$basename.org_counts') );
my $t4 = time();
printf STDERR "... %.3f seconds\n\n", $t4-$t3;


