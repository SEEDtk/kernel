#! /usr/bin/env perl

use strict;

use Carp;
use Data::Dumper;
use Getopt::Long;
use gjoseqlib;

my ($help, $len, $cov, $mean, $stdev, $err, $all2all);

GetOptions("h|help" => \$help,
           "a"      => \$all2all,
           "c=f"    => \$cov,
           "e=f"    => \$err,
           "l=i"    => \$len,
           "m=f"    => \$mean,
           "s=f"    => \$stdev,
	  ) or die("Error in command line arguments\n");

my $usage = <<End_Usage;

Usage: $0 [options] ref.fa

       -h          - print this help statement
       -a          - show all to all matches in reference vs assembly dot plot
       -c cov      - fold coverage used in simulating reads (D: 20)
       -e fract    - fraction of standard Illumina HiSeq 2000 error rate (D: 1.0)
       -l len      - simulated read length (D: 100)
       -m mean     - mean of paired end read fragment size (D: 300)
       -s stdev    - standard deviation of read fragment size (D: 10)

End_Usage

$help and die $usage;

my $ref = shift @ARGV or die $usage;

my $reads = sim_reads($ref, $len, $cov, $mean, $stdev, $err);

assemble_with_velvet($reads, $ref, $all2all);
assemble_with_spades($reads, $ref, $all2all);

print STDERR "Done.\n";

sub assemble_with_spades {
    my ($reads, $ref, $all2all) = @_;
    my $out     = output_prefix($ref);
    my $spades  = "/home/fangfang/bin/spades.py";
    my $tmpdir  = "tmp.$out.spades";
    my $contigs = "$out.spades.fa";

    print STDERR "Assembling reads with spades into '$contigs'...\n";
    run("$spades -1 $reads->[0] -2 $reads->[1] -o $tmpdir &>spades.log");
    run("cp $tmpdir/contigs.fasta $out.spades.raw.fa");
    run("rm -rf $tmpdir");
    filter_contigs("$out.spades.raw.fa", $contigs);
    compare_ref_and_assembly($ref, $contigs, $all2all);
}

sub assemble_with_velvet {
    my ($reads, $ref, $all2all) = @_;
    my $out     = output_prefix($ref);
    my $velvetg = "/home/fangfang/bin/velvetg";
    my $velveth = "/home/fangfang/bin/velveth";
    my $tmpdir  = "tmp.$out.velvet";
    my $contigs = "$out.velvet.fa";

    print STDERR "Assembling reads with velvet into '$contigs'...\n";
    run("$velveth $tmpdir 23 -shortPaired -fastq -separate $reads->[0] $reads->[1] &>velvet.log");
    run("$velvetg $tmpdir &>velvet.log");
    run("cp $tmpdir/contigs.fa $out.velvet.raw.fa");
    run("rm -rf $tmpdir");
    filter_contigs("$out.velvet.raw.fa", $contigs);
    compare_ref_and_assembly($ref, $contigs, $all2all);
}

sub compare_ref_and_assembly {
    my ($ref, $assembly, $all2all) = @_;
    my $out    = output_prefix($assembly);
    my $nucmer = "/home/fangfang/bin/nucmer";
    my $plot   = "/home/fangfang/bin/mummerplot";
    my $webdir = "$ENV{HOME}/public_html/CovCov";
    my $url    = "http://bioseed.mcs.anl.gov/~$ENV{USER}/CovCov/$out.png";
    my $filter = $all2all ? '' : '--filter';
    run("$nucmer -p $out $ref $assembly --maxmatch --coords &>mummer.log");
    run("$plot --title 'REF vs Assembly ($ref vs $assembly)' -p $out --png $out.delta $filter --layout -R $ref -Q $assembly &>mummer.log");
    run("mkdir -p $webdir");
    run("cp $out.png $webdir");
    print STDERR "Plotted dot plot: $url\n";
    -e "$out.$_" && unlink "$out.$_" for qw(fplot rplot filter gp);
}

sub filter_contigs {
    my ($in, $out, $minlen, $mincov) = @_;

    $minlen ||= 150;
    $mincov ||= 5;

    my @seqs = gjoseqlib::read_fasta($in);
    @seqs = grep { length $_->[2] >= $minlen } @seqs;

    if ($mincov > 0) {
        @seqs = grep { ($_->[0].$_->[1]) !~ /cov[_ :=]([0-9.]+)/i || $1 >= $mincov } @seqs;
    }

    gjoseqlib::print_alignment_as_fasta($out, \@seqs);
}

sub sim_reads {
    my ($ref, $len, $cov, $mean, $stdev, $err) = @_;

    $len   ||= 100;
    $cov   ||= 20;
    $mean  ||= 300;
    $stdev ||= 10;
    $err     = 1.0  if !defined($err);
    $err     = 1e-8 if $err < 1e-8;

    print STDERR "Generating simulated reads with art_illumina...\n";
    print STDERR "    read length   = $len\n";
    print STDERR "    fold coverage = $cov\n";
    print STDERR "    insert size   = $mean (±$stdev) bp\n";
    print STDERR "    error rate    = ${err}X Illumina HS2000 error (~1%)\n";

    my $qs = -10 * log($err)/log(10);
    my $art = "/home/fangfang/bin/art_illumina";
    my $out = output_prefix($ref).'-read';
    my $cmd = "$art -ss HS20 -i $ref -l $len -f $cov -p -m $mean -s $stdev -qs $qs -qs2 $qs -na -o $out &>art.log";

    run($cmd);
    my @reads = map { "$out$_.fq" } (1, 2);

    wantarray ? @reads : \@reads;
}

sub output_prefix {
    my ($file) = @_;
    $file =~ s/^.*\///;
    $file =~ s/\.(fa|fasta)$//;
    $file;
}

sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
