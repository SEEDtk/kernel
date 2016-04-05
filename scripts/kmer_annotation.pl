#!/usr/bin/perl

use strict;
use Data::Dumper;

my $kmerH = {};
my ($training,$test,$projections);

my $k = 8;

my $usage = "usage: perl kmer_annotation.pl TrainingF TestF kmer-size > ProjectionF\n";
(
 ($training =    shift @ARGV) &&
 ($test     =    shift @ARGV) &&
 ($k        =    shift @ARGV)
)
    || die $usage;

open(TRAINING,"<$training") || die "could not open $training";
while (defined($_ = <TRAINING>))
{
    chomp;
    my($peg,$func,$seq) = split(/\t/,$_);
    my $i;
    for ($i=0; ($i < (length($seq) - $k)); $i++)
    {
	my $kmer = uc substr($seq,$i,$k);
        $kmerH->{$kmer}->{$func}++;
    }
}
close(TRAINING);
open(TEST,"<$test") || die "could not open $test";
while (defined($_ = <TEST>))
{
    chomp;
    my($peg,$func,$seq) = split(/\t/,$_);

    &process1($peg,$func,$seq,$kmerH);
#   &process2($peg,$func,$seq,$kmerH);
}

sub process1 {
    my($peg,$func,$seq,$kmerH) = @_;

    my %counts;
    my $i;
    for ($i=0; ($i < (length($seq) - $k)); $i++)
    {
	my $kmer = uc substr($seq,$i,$k);
	my $subH = $kmerH->{$kmer};
        my @keys = keys(%$subH);
	if (@keys == 1)
	{
	    $counts{$keys[0]}++;
	}
    }
    my @poss = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
    my $f;
    if ((@poss == 1) || (@poss > 1) && ($counts{$poss[0]} >= (2 * $counts{$poss[1]})))
    {
	$f = $poss[0];
    }
    else
    {
	$f = "hypothetical protein";
    }
    print join("\t",($peg,$func,$seq,$f)),"\n";
}

sub process2 {
    my($peg,$func,$seq,$kmerH) = @_;

    my %counts;
    my $i;
    for ($i=0; ($i < (length($seq) - $k)); $i++)
    {
	my $kmer = uc substr($seq,$i,$k);
	my $subH = $kmerH->{$kmer};
        my @funcs = keys(%$subH);
        foreach my $f (@funcs)
	{
	    $counts{$f} += $subH->{$f};
	}
    }
    my @poss = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
    my $f;
    if ((@poss == 1) || (@poss > 1) && ($counts{$poss[0]} >= (2 * $counts{$poss[1]})))
    {
	$f = $poss[0];
    }
    else
    {
	$f = "hypothetical protein";
    }
    print join("\t",($peg,$func,$seq,$f)),"\n";
}
