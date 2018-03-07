#!/usr/bin/env perl

=head1 Create Probdir: Function Predictors Step 2

    build_matrix raw.table probDir roles.to.use

This script takes the output from L<build_role_tables.pl> and creates a function matrix. Each row of the
matrix represents a genome. Each column of the matrix represents a role. The values in the matrix indicate
how many times the role occurs in the genome. The resulting matrix is encoded as a SciKit probdir. In other
words, there is a C<row.h> of genome IDs, an C<col.h> of role IDs, and an C<X> file containing the actual matrix.

This is the second step of the process for building function predictors. The third is L<build_predictor_set.pl>.

=head2 Parameters

The positional parameters are the name of the C<raw.table> file, the name of the output directory, and a
filter file. If specified, the filter file contains a list of the permissible roles. (If used, this should be
a subset of the C<roles.in.subsystems> output by L<build_role_tables.pl>).

The output directory cannot already exist unless the C<--clear> option is specified.

=cut

use strict;
use warnings;
use Data::Dumper;
use ScriptUtils;
use SeedUtils;
use File::Copy::Recursive;

my %funcs;
my %genomes;
my %counts;

my $opt = ScriptUtils::Opts('raw.table probDir roles.to.use',
        ['clear', 'overwrite previous results']);

my ($table_file, $probDir, $keep) = @ARGV;

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
} elsif ($opt->clear) {
    print "Clearing $probDir.\n";
    if (-d "$probDir/Predictors") {
        File::Copy::Recursive::pathrmdir("$probDir/Predictors") || die "Could not clear Predictors directory: $!";
    }
    my @files = qw(col.h roles.mapped roles.not_mapped row.h X);
    for my $file (@files) {
        if (-f "$probDir/$file") {
            unlink "$probDir/$file";
        }
    }
} else {
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
    my (undef, $funcID, $fid) = split /\t/, $line;

    if ($keep) {
        next unless $funcs{$funcID};
    }
    else {
        $funcs{$funcID} = 1;
    }

    my $genome = &SeedUtils::genome_of($fid);
    $genomes{$genome} = 1;

    ++$counts{$genome}->{$funcID} ;
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
