#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use SeedUtils;

my %funcs;
my %genomes;
my %counts;

my ($table_file, $probDir, $keep) = @ARGV;

my %keep;
if ($keep) {
    if (-s $keep) {
	%funcs = map { m/^(\S+)/ ? ($1 => 1) : () } &SeedUtils::file_read($keep);
    }
    else {
	die "Filter-file '$keep' does not exist";
    }
}

if (!-d $probDir) {
    mkdir($probDir) or die "Could not create probDir '$probDir'";
}
else {
    die "ERROR: probDir '$probDir' already exists";
}

my $table_fh;
open($table_fh, '<', "$table_file")    or die "Could not read-open '$table_file'";

my ($X_fh, $row_fh, $col_fh);
open($X_fh,     '>', "$probDir/X")     or die "Could not write-open '$probDir/X'";
open($row_fh,   '>', "$probDir/row.h") or die "Could not write-open '$probDir/row.h'";
open($col_fh,   '>', "$probDir/col.h") or die "Could not write-open '$probDir/col.h'";

my $line;
my $tick = 0;
my $tock = 0;
while (defined($line = <$table_fh>)) {
    if ($tick >= 10000) { print STDERR ".";  $tick = 0; $tock++; } else { $tick++; }
    if ($tock >= 100)   { print STDERR "\n"; $tock = 0; }
    
    chomp $line;
    my (undef, $func, $fid) = split /\t/, $line;
    
    if ($keep) {
	next unless $funcs{$func};
    }
    else {
	$funcs{$func} = 1;
    }
    
    my $genome = &SeedUtils::genome_of($fid);
    $genomes{$genome} = 1;
    
    ++$counts{$genome}->{$func} ;
}
print STDERR "\n\n";

my @genomes = sort keys %genomes;
my @funcs   = sort keys %funcs;

for (my $i=0; $i < @genomes; ++$i) {
    print $row_fh ($i, "\t", $genomes[$i], "\n");
}

for (my $j=0; $j < @funcs; ++$j) {
    print $col_fh ($j, "\t", $funcs[$j], "\n");
}

for (my $i=0; $i < @genomes; ++$i) {
    print $X_fh (join("\t", map { sprintf('%.1f', $counts{$genomes[$i]}->{$_} || 0.0) } @funcs), "\n");
}
