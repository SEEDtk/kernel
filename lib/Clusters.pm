package Clusters;

use strict;
use Data::Dumper;
use BasicLocation;

sub subsys_mems {
    my($genome,$shrub) = @_;
    my @tuples = $shrub->GetAll("Genome Genome2Row SubsystemRow Row2Cell SubsystemCell Cell2Feature",
                                "Genome(id) = ?",
                                [$genome],
                                "Cell2Feature(to-link) SubsystemRow(id) SubsystemRow(variant-code)");
    my $peg_to_subsys_row = {};
    foreach my $tuple (@tuples)
    {
        my($peg,$ss_row,$vc) = @$tuple;
        if (&SeedUtils::probably_active($vc))
        {
            $peg_to_subsys_row->{$peg}->{$ss_row} = 1;
        }
    }
    return $peg_to_subsys_row;
}

sub cluster_subsys_fids {
    my($max_gap,$sorted_fid_tuples,$peg_to_subsys_row) = @_;

    my %pairs;
    my $i=0;
    while ($i < (@$sorted_fid_tuples - 1))
    {
        my @ss = &clustered_pair($sorted_fid_tuples->[$i],$sorted_fid_tuples->[$i+1],$max_gap,$peg_to_subsys_row);
        if (@ss > 0)
        {
            my $fid1 = $sorted_fid_tuples->[$i]->[0];
            my $fid2 = $sorted_fid_tuples->[$i+1]->[0];
            foreach my $ss (@ss)
            {
                push(@{$pairs{$ss}},[$fid1,$fid2]);
            }
        }
        $i++;
    }

    my @ss_clusters;
    foreach my $ss (keys(%pairs))
    {
        my $con = $pairs{$ss};
        my @clusters = &cluster_connections($con);
        foreach my $cluster (@clusters)
        {
            push(@ss_clusters,[$ss,$cluster]);
        }
    }
    return @ss_clusters;
}

sub clustered_pair {
    my($tuple1,$tuple2,$max_gap,$peg_to_subsys_row) = @_;

    my @ss = ();
    my ($fid1,$loc1) = @$tuple1;
    my ($fid2,$loc2) = @$tuple2;
    if ($loc1->Contig ne $loc2->Contig) { return () }
    my $midpt1 = ($loc1->Left + $loc1->Right) / 2;
    my $midpt2 = ($loc2->Left + $loc2->Right) / 2;
    if (abs($midpt1 - $midpt2) > $max_gap) { return () }
    my $ss1 = $peg_to_subsys_row->{$fid1};
    my $ss2 = $peg_to_subsys_row->{$fid2};
    if ($ss1 && $ss2)
    {
        @ss = &intersection($ss1,$ss2);
    }
    return @ss;
}

sub intersection {
    my($xH,$yH) = @_;
    return grep { $yH->{$_} } keys(%$xH);
}

sub cluster_connections {
    my($pairs) = @_;
    my %to_cluster;
    my %in_cluster;
    my $nxt = 1;

    foreach my $pair (@$pairs)
    {
        my($obj1,$obj2) = @$pair;
        my $in1 = $to_cluster{$obj1};
        my $in2 = $to_cluster{$obj2};

        if (defined($in1) && defined($in2) && ($in1 != $in2))
        {
            push(@{$in_cluster{$in1}},@{$in_cluster{$in2}});
            foreach $_ (@{$in_cluster{$in2}})
            {
                $to_cluster{$_} = $in1;
            }
            delete $in_cluster{$in2};
        }
        elsif ((! defined($in1)) && defined($in2))
        {
            push(@{$in_cluster{$in1}},@{$in_cluster{$in2}});
            foreach $_ (@{$in_cluster{$in2}})
            {
                $to_cluster{$_} = $in1;
            }
            delete $in_cluster{$in2};
        }
        elsif ((! defined($in1)) && defined($in2))
        {
            push(@{$in_cluster{$in2}},$obj1);
            $to_cluster{$obj1} = $in2;
        }
        elsif ((! defined($in2)) && defined($in1))
        {
            push(@{$in_cluster{$in1}},$obj2);
            $to_cluster{$obj2} = $in1;
        }
        elsif ((! defined($in1)) && (! defined($in2)))
        {
            $to_cluster{$obj1} = $to_cluster{$obj2} = $nxt;
            $in_cluster{$nxt} = [$obj1,$obj2];
            $nxt++;
        }
    }
    return map { $in_cluster{$_} } keys(%in_cluster);
}

1;
