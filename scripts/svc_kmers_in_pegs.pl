#!/usr/bin/perl

use strict;
use Data::Dumper;
use Shrub;
use ScriptUtils;
use File::Copy::Recursive;

=head1 Generate Classification Directory for Kmers in Pegs

    svc_kmers_in_pegs -d ClassDir < pegs

Takes a file of pegs as input.  Creates a directory

    ClassDir
        X
        y
        col.h
        row.h
        
We are creating a matrix that indicates presence or absence of kmers in pegs. Each
row is a peg. Each column is a kmer.

=head2 Parameters

=over 4

=item ClassDir

The name of a generated directory used to classify pegs based on Kmer content

=back

=cut

$| = 1;
my $opt = ScriptUtils::Opts( '',
                 Shrub::script_options(), 
                 ScriptUtils::ih_options(),
                 ['col|c=i', 'pegid column', { }],
                 ['ksize|k=i','kmer size',{ default => 6 } ],
                 ['dir|d=s', 'Name of generated classification directory',{required => 1}],
                 ['force|f', 'Delete old copy of work directory.']
        );

# Connect to the database.
my $shrub  = Shrub->new_for_script($opt);
my $ih     = ScriptUtils::IH( $opt->input );
my $column = $opt->col;
my $dir    = $opt->dir;
my $k      = $opt->ksize;

if (-d $dir) {
    if ($opt->force) {
        File::Copy::Recursive::pathempty($dir) || die "Could not empty previous contents of $dir.";
    } else {
        die "$dir already exists.";
    }
} else {
    mkdir($dir, 0777);
}

# X is the main array. One row for each peg, containing 1s and 0s that indicate presence or absence of a kmer
# col.h are the column headrs (indexes of kmers)
# row.h are the row headers (indexes of pegs)
# pred.h are the function index (indexes of functions)
# 
# y is the array of dependent variables (the predictions), mapping pegs to function indices. This file contains
#   one number per row. The row corresponds to the peg in row.h. The number is the index of the function in pred.h

open(X,">$dir/X") || die "could not open $dir/X";
open(Y,">$dir/y") || die "could not open $dir/y";
open(COL,">$dir/col.h") || die "could not open $dir/col.h";
open(ROW,">$dir/row.h") || die "could not open $dir/row.h";
open(FUNCS,">$dir/pred.h") || die "could not open $dir/pred.h";

my $kmerN = 0;
my $funcN = 0;
my $pegN  = 0;
# This hash maps a kmer to its index in the col.h file.
my $kmerH = {};
# This is a two-level hash by peg ID and kmer that indicates which kmers occur in each peg.
my %kmer_occurs;
# This hash maps each function to its index in the pred.h file.
my $funcI = {};
# This is the final list of kmers.
my @kmers;
# This is the final list of pegs.
my @pegs;
# Loop through the batches.
my $pegtot = 0;
while (my @tuples = ScriptUtils::get_couplets($ih, $column, 50000)) {
    my $pegNow = scalar(@tuples);
    $pegtot += $pegNow;
    print "$pegNow read in this batch. $pegtot total.\n";
    my @ids = map { $_->[0] } @tuples;
    my $seqH = $shrub->Feature2Trans(\@ids);
    print scalar(keys %$seqH) . " translations found.\n";
    my $funcH = $shrub->Feature2Function(1,\@ids);
    print scalar(keys %$funcH) . " functions found.\n";


    foreach my $id (keys(%$seqH))
    {
        print ROW $pegN++,"\t$id\n";
        push @pegs, $id;
        my $seq = $seqH->{$id};
        my $func = $funcH->{$id}->[1];
        if (! defined $funcI->{$func})
        {
            print FUNCS $funcN,"\t",$func,"\n";
            $funcI->{$func} = $funcN++;
        }
        print Y $funcI->{$func},"\n";
        my $i;
        for ($i=0; ($i < (length($seq) - $k)); $i++)
        {
            my $kmer = uc substr($seq,$i,$k);
            if (! $kmerH->{$kmer})
            {
                push(@kmers,$kmer);
                $kmerH->{$kmer} = $kmerN++;
            }
            $kmer_occurs{$id}->{$kmer} = 1;
        }
    }
}
close(ROW);
close(Y);
close(FUNCS);
print "spooling kmers.\n";
for (my $i=0; ($i < @kmers); $i++)
{
    print COL join("\t",($i,$kmers[$i])),"\n";
}
close(COL);
my $pegCount = 0;
foreach my $id (@pegs)
{
    my @row;
    for (my $j = 0; ($j < @kmers); $j++)
    {
        my $v = $kmer_occurs{$id}->{$kmers[$j]};
        push(@row,$v ? 1 : 0);
    }
    print X join(",",@row),"\n";
    $pegCount++;
    if ($pegCount % 50000 == 0) {
        print "$pegCount pegs output.\n";
    }
}

close(X);
print "All done. $pegCount total pegs.\n";