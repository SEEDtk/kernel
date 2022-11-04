#!/usr/bin/env perl
#
#  overlap_reads  fwd_read_file  rev_read_file  > merged_reads
#
use strict;
use gjoseqlib;
use Data::Dumper;

use P3Utils;

=head1 Merge left- and right-end reads that sufficiently overlap

    overlap_reads  [options]  fwd_read_file  rev_read_file  > merged_reads

This script merges left- and right-end reads into a single spanning equence
if they overlap by more than a specified amount and are sufficiently identitcal
in the overlap region

=head2 Parameters

There are two positional parameters:

=over 4

=item fwd_read_file

FASTA file of "forward"  (AKA "left-end") reads.


=item rev_read_file

FASTA file of "reverse"  (AKA "right-end") reads.

=back

The following command-line parameters are supported:

=over 4

=item probe (p)

The minimum length of contiguously matching characters within the overlap region.
(Default: 25nt)

=item strict (s)

Require 100% identity over the full length of the overlap region.
(Default: "non-strict" --- require 100% identity only over the width of the probe-length)

=item verbose

If specified print detailed information to error filehandle.

=item debug                                                                                                                    

If specified print debugging information to error filehandle.

=item help (h)

Print usage-message and exit.

=back

=cut 

my $usage = <<'End_of_usage';

Usage: overlap_reads  [options]  fwd_read_file  rev_read_file  > merged_reads

Options:

  -h -help        #  Print usage information and quit
  -p -probe  int  #  The required length of perfect matching (D = 25) 
  -s -strict      #  Demand identity in overlap of forward and reverse reads;
                  #      default is identity over the probe length defined by -p

End_of_usage


my $opt = P3Utils::script_opts(
    'fwd_read_file  rev_read_file',
    [ 'probe|p=i',  'The minimum required length for a "perfect match" ovwrlap-region (D: 25)',   { default => 25 } ],
    [ 'strict|s',   'Demand identity over the full overlap-region betwwee the forward and reverse reads. (D: "non-strict" --- require identity only over the probe length)',
      { default => 0 } ],
    [ 'verbose',    'Print detailed information to error-filehandle',  ],
    [ 'debug',      'Print debugging information to error-filehandle', ],
    );

my $strict    = $opt->strict();    #  Default: 0   (Do not demand perfect match in overlap)
my $probe_len = $opt->probe();     #  Default: Arbitrarily demand a perfect overlap of 20 nt

my $verbose   = $opt->verbose() || $ENV{VERBOSE};
my $debug     = $opt->debug()   || $ENV{DEBUG};


print STDERR Dumper(\@ARGV) if $debug;
my ( $fwd_file, $rev_file ) = @ARGV;
print STDERR ("fwd_file=$fwd_file\n", "rev_file=$rev_file\n") if $debug;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Check if mandatory input-files exist, and open for read-access...
#-----------------------------------------------------------------------
-s $fwd_file && -s $rev_file
    or print STDERR "Could not find sequence files.\n", $usage
        and exit;

open( FWD, '<', $fwd_file )
    or die "Could no open forward read file '$fwd_file'";
open( REV, '<', $rev_file )
    or die "Could no open reverse read file '$rev_file'";

my ( $fwd_header, $rev_header );

my $pairs  = 0;  #  Read pairs processed
my $merge1 = 0;  #  Merged in first try
my $merge2 = 0;  #  Merged in second try
my @length;

while ( defined( $fwd_header = <FWD> )
     && defined( $rev_header = <REV> )  #  Not used
      )
{
    my $fwd_seq = <FWD>;
    my $rev_seq = <REV>;
    last unless $fwd_seq && $rev_seq;

    chomp $fwd_seq;
    chomp $rev_seq;
    $rev_seq = reverse $rev_seq;
    $rev_seq =~ tr/acgtACGT/tgcaTGCA/;

    $pairs++;

    #  Prefix of reverse read matches forward read
    #
    #                           offset1
    #                              |
    #   forward   ----------------------------------
    #   reverse                    |-----|------------------------------
    #                               probe
    #
    #  Suffix of forward read matches reverse read
    #
    #                                         probe
    #   forward   ---------------------------|-----|
    #   reverse                    -------------------------------------
    #                                        |
    #                                     offset2
    #
    my $offset1;
    my $offset2;
    my $merged;
    if ( $strict )
    {
        my $offset1 = index( $fwd_seq, substr( $rev_seq, 0, $probe_len ) );
        if ( $offset1 >= 0 )
        {
            my $overlap1 = substr( $fwd_seq, $offset1 );
            my $overlap2 = substr( $rev_seq, 0, length( $overlap1 ) );
            if ( $overlap1 eq $overlap2 )
            {
                $merged = substr( $fwd_seq, 0, $offset1 ) . $rev_seq;
                $merge1++;
            }
        }
    }
    else
    {
        my $offset1 = index( $fwd_seq, substr( $rev_seq, 0, $probe_len ) );
        if ( $offset1 > 0 )
        {
            $merged = substr( $fwd_seq, 0, $offset1 ) . $rev_seq;
            $merge1++;
        }

        #  This is a fallback if the end of the reverse read does not match in forward read.
        #  We should probably throw out all imperfect overlaps.

        elsif ( ( $offset2 = index( $rev_seq, substr( $fwd_seq, -$probe_len ) ) ) > 1 )
        {
            $merged = $fwd_seq . substr( $rev_seq, $offset2+$probe_len );
            $merge2++;
        }
    }

    print "$fwd_header$merged\n"  if $merged;  #  $fwd_header still has its new line
    $length[length($merged)]++    if $merged;
}

printf STDERR "%9d pairs read\n",           $pairs;
printf STDERR "%9d merged\n",               $merge1+$merge2;
if ( ! $strict )
{
    printf STDERR "%9d merged in first try\n",  $merge1;
    printf STDERR "%9d merged in second try\n", $merge2;
}
printf STDERR "%9d failed\n",               $pairs-($merge1+$merge2);

print STDERR "\n";
print STDERR "Length distribution of merged reads:\n";
print STDERR "---------------\n";
print STDERR "Length    Count\n";
print STDERR "---------------\n";

my $print;
for ( my $len = 0; $len < @length; $len++ )
{
    my $cnt = $length[$len] || 0;
    $print ||= $cnt;
    printf STDERR "%4d%9d\n", $len, $cnt  if ($print && ($cnt || $debug));
}

exit(0);
