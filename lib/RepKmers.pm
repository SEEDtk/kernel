package RepKmers;

use strict;
use Data::Dumper;

sub sim {
    my($s1,$s2,$K) = @_;

    my %in1;
    my $common = 0;
    my $i;
    for ($i=0; ($i < (length($s1)-$K)); $i++)
    {
	$in1{lc substr($s1,$i,$K)} = 1;
    }

    for ($i=0; ($i < (length($s2)-$K)); $i++)
    {
	if ($in1{lc substr($s2,$i,$K)})
	{
	    $common++;
	}
    }
    return $common;
}

sub kmers_of_seq {
    my($K,$seq) = @_;

    my $kmers = {};
    my $i;
    for ($i=0; ($i < (length($seq)-$K)); $i++)
    {
	$kmers->{lc substr($seq,$i,$K)} = 1;
    }
    return $kmers;
}

sub in_common {
    my($K,$kmers,$seq) = @_;

    my $common = 0;
    my $i;
    for ($i=0; ($i < (length($seq)-$K)); $i++)
    {
	if ($kmers->{lc substr($seq,$i,$K)})
	{
	    $common ++;
	}
    }
    return $common;
}

1;
