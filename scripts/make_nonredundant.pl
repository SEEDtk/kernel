#!/usr/bin/env perl
#
#  make_nonredundant < merged_reads > unique_reads
#
#  Note:  There is now a default setting that eliminates sequences that are
#         less than 0.0005% of the total reads.  See the -f and -n options
#         to 
#
use strict;
use gjoseqlib;
use Data::Dumper;

use P3Utils;

# my $usage = <<'End_of_usage';
#
# Usage: make_nonredundant [options] < merged_reads > unique_reads
#
# Options:
#
#   -c -count     int    #  Same is -min_count
#   -f -min_freq  fract  #  Minimum frequency of occurance to retain a sequence
#                        #      (D = 5e-6, i.e., 0.0005%)
#   -h -help             #  Print usage information and quit
#   -n -min_count int    #  The minimum number of occurances to keep a sequence
#                        #      (by default, 5e-6 * number of reads)
#
# The min_count option takes presidence over min_freq.
#
# End_of_usage

=head1 Removes duplicate sequences from a FASTA file.

    make_nonredundant [options] < merged_reads > unique_reads

This script scans an input FASTA file, counts the number of instances of each sequence,
and outputs a single instance of each sequence that exceeds a minimum number or minimum frequency.

=head2 Parameters                                                                                                                                                                                                                                            
There are no positional parameters.                                                                                            

The following command-line parameters are supported:        

=over 4

=item count (c)

Minimum number of instances to keep a sequence. (Default: 5e-6 * number of reads)

=item min-count (n)

Minimum number of instances to keep a sequence. (Default: 5e-6 * number of reads)

item min-freq (f)

Minimum frequency of occurance to keep a sequence.  (Deafault: 5e-6)

=item help (h)

Display usage information and exit.

=back

=cut

my ($opt, $usage) = P3Utils::script_opts(
    '', P3Utils::ih_options(),
    [ 'output|o=s', 'Output file (if not STDOUT)', ],
    [],
    [ 'count|c=i',     'Minimum number of instances to keep a sequence. (D: 5e-6 * number of reads)',
      { default => undef } ],
    [ 'min-count|n=i', 'synonym of count|c',
      { default => undef } ],
    [ 'min-freq|f=i',  'Minimum frequency of occurance to keep a sequence.  (D: 5e-6)',
      { default => 5.0e-5 } ],
    );

my $min_cnt   = $opt->min_count() || $opt->count();   # 0.0005%
my $min_freq  = $opt->min_freq();                     # 0.0005%

my $help        = $opt->help();
my $input_file  = $opt->input();
my $output_file = $opt->output();


my $input_fh  = \*STDIN;
if ($input_file) {
    open($input_fh, q(<), $input_file)
	or die("Could not read-open input-file '$input_file'");
}

my $output_fh = \*STDOUT;
if ($output_file) {
    open($output_fh, q(>), $output_file)
	or die("Could not write-open output-file '$output_file'");
}


# GetOptions( "c|count=i"     => \$min_cnt,
#             "f|min_freq=f"  => \$min_freq,
#             "h|help"        => \$help,
#             "n|min_count=i" => \$min_cnt
#           )
#     or print STDERR "Error in command line arguments.\n", $usage
#         and exit;

if ( $help ) {
    print STDERR $usage;
    exit;
}

my %cnt;
my $n;
while (<$input_fh>)
{
    next if /^>/;
    chomp;
    $cnt{$_}++;
    $n++;
}

$min_cnt = int( $n * $min_freq )  if ! defined( $min_cnt );
my $seq_omit  = 0;
my $read_omit = 0;

my $n_seq = keys(%cnt);
my @data = sort { $b->[1] <=> $a->[1] || $a->[0] cmp $b->[0] }
           map  { my $cnt = $cnt{$_};
                  if ( $cnt < $min_cnt )
                  {
                      $read_omit += $cnt;
                      $seq_omit++;
                  }
                  $cnt >= $min_cnt ? [ $_, $cnt ] : ();
               }
           keys %cnt;

my $seq_kept  = 0;
my $read_kept = 0;
foreach ( @data )
{
    my $cnt = $_->[1];
    $read_kept += $cnt;
    printf $output_fh ">seq%08d count = %d\n", ++$seq_kept, $cnt;
    printf $output_fh "$_->[0]\n";
}

print STDERR  "$n_seq unique sequences found among $n merged reads\n\n";
if ( $seq_omit )
{
    printf STDERR "Sequences were omitted due to low number of occurrences (< %d reads, or %.4f%% of the total):\n", $min_cnt, 100*$min_cnt/$n;
    print  STDERR "You can change this behavior with the -f or -n option.\n";
    printf STDERR "Omitted:  %8d sequences (%7.4f%%), representing %8d reads (%7.4f%%).\n", $seq_omit, 100*$seq_omit/$n_seq, $read_omit, 100*$read_omit/$n;
    printf STDERR "Retained: %8d sequences (%7.4f%%), representing %8d reads (%7.4f%%).\n", $seq_kept, 100*$seq_kept/$n_seq, $read_kept, 100*$read_kept/$n;
    print  STDERR "\n";
}

