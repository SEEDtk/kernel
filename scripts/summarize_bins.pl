use strict;
use Data::Dumper;

$/ = "\n//\n";

while (defined ($_ = <STDIN>))
{
    $/ = "\n";
    # Process the data for a bin.
    &process_1($_);
    $/ = "\n//\n";
}

sub process_1 {
    my($x) = @_;
    # Split the universal roles and the contig list.
    if ($x =~ /^(.*)\n\n(.*)\n?\/\/\n/s)
    {
        my $univs = $2;
        print $univs,"\n--------\n";
        my $contigs = $1;
        my @counts = &count_ref($contigs);
        foreach my $tuple (@counts)
        {
            print join("\t",@$tuple),"\n";
        }
        print "\n--------\n";
        my $sz = &sum_lengths($contigs);
        print "total size of contigs=$sz\n";
        print "//\n\n";
    }
}


sub count_ref {
    my($contigs) = @_;

    my %counts;
    my %bp;
    my @entries = split(/\n\n/,$contigs);
    foreach my $x (@entries)
    {
        # Get the contig length and the genome name.
        if ($x =~ /length_(\d+).*\t([^\t]*)$/s)
        {
            $bp{$2} += $1;
            $counts{$2}++;
        }
    }

    return sort { $b->[0] <=> $a->[0] } map { [$counts{$_},$bp{$_},$_] } keys(%counts);
}

sub sum_lengths {
    my($contigs) = @_;

    my $tot=0;
    my @entries = split(/\n\n/,$contigs);
    foreach my $x (@entries)
    {
        if ($x =~ /length_(\d+)/)
        {
            $tot += $1;
        }
    }
    return $tot;
}
