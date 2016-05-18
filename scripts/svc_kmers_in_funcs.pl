#!/usr/bin/perl

use strict;
use Data::Dumper;
use Shrub;
use ScriptUtils;
use File::Copy::Recursive;

=head1 Generate Classification Directory for Kmers in Functions

    svc_kmers_in_funcs -d ClassDir < pegs

Takes a file of pegs as input.  Creates a directory

    ClassDir
        X
        y
        col.h
        row.h
        y.map

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
                             ['force|f', 'Delete old copy of work directory.'],
			     ['dir|d=s', 'Name of generated classification directory',{required => 1}]
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

# X is the main array.
# col.h are the column headrs (indexes of kmers)
# row.h are the row headers (indexes of pegs)
# y.map are the function index
# 
# y is the array of dependent variables (the predictions)

open(X,">$dir/X") || die "could not open $dir/X";
open(Y,">$dir/y") || die "could not open $dir/y";
open(COL,">$dir/col.h") || die "could not open $dir/col.h";
open(ROW,">$dir/row.h") || die "could not open $dir/row.h";

my $kmerN = 0;
my $funcN = 0;
my %kI;
my %fI;
my $i;
my %kmer_occurs;
my $count = 0;
my @pegs;
while (defined($_ = <$ih>))
{
    chomp;
    my @fields = split(/\t/,$_);
    if (! $column)
    {
	push(@pegs,$fields[-1]);
    }
    else 
    {
	push(@pegs,$fields[$column-1]);
    }
    $count++;
    if ($count % 50000 == 0) {
        print "$count pegs processed.\n";
    }
}
print "$count pegs total.\n";
my $seqH  = $shrub->Feature2Trans(\@pegs);
my $funcH = $shrub->Feature2Function(1,\@pegs);
my %funcs = map { $_->[1] => 1 } grep { defined $_ } values %$funcH;
my @functions = sort keys %funcs;               # sorted list of functions; rows of X
%funcs = ();
$count = 0;
for ($i=0; ($i < @functions);$i++)
{
    $fI{$functions[$i]} = $i;                   # $fI maps functions to rows
    print ROW $i,"\t",$functions[$i],"\n";
    print Y $i,"\n";                            # one Y value per function
    $count++;
    if ($count % 50000 == 0) {
        print "$count functions processed.\n";
    }
}
close(ROW);
close(Y);
print "$count total functions.\n";
my @kmers;
$count = 0;
foreach my $peg (keys(%$seqH))
{
    my $seq = $seqH->{$peg};
    if ($seq) {
        my $func = $funcH->{$peg}->[1];
        my $funcI = $fI{$func};                      # funcI is the row corresponding to the function
        my $i;
        for ($i=0; ($i < (length($seq) - $k)); $i++)
        {
            my $kmer = uc substr($seq,$i,$k);
            if (! defined($kI{$kmer}))
            {
                push(@kmers,$kmer);
                $kI{$kmer} = $kmerN++;
            }
            $kmer_occurs{$funcI}{$kI{$kmer}} = 1;
        }
        $count++;
        if ($count % 50000 == 0) {
            print "$count sequences processed.\n";
        }
    }
}
$count = 0;
for (my $i=0; ($i < @kmers); $i++)
{
    print COL join("\t",($i,$kmers[$i])),"\n";
    $count++;
    if ($count % 50000 == 0) {
        print "$count kmers processed.\n";
    }
}
close(COL);
print "$count total kmers.\n";
my $zeros = join(",", map { '0' } @kmers);
$count = 0;
for ($i=0; ($i < @functions); $i++)
{
    my $row = $zeros;
    my $occurs = $kmer_occurs{$i};
    for my $j (keys %$occurs) {
        substr $row, 2*$j, 1, '1';
    }
    print X "$row\n";
    $count++;
    if ($count % 50000 == 0) {
        print "$count rows processed.\n";
    }
}
close(X);
print "Copying row.h to y.map\n";
File::Copy::Recursive::fcopy("$dir/row.h", "$dir/y.map");
print "All done. $count total rows produced.\n";
