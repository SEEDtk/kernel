package Triples;

use strict;
use Data::Dumper;
use BasicLocation;

sub subsystem_roles {
    my($shrub) = @_;

    my @tuples = $shrub->GetAll("Subsystem2Role Role", "", [],
				"Subsystem2Role(to-link) Role(description)");
    my $roles = {};
    foreach my $tuple (@tuples)
    {
	my($role,$desc) = @$tuple;
	$roles->{$role} = $desc;
    }
    return $roles;
}

sub triples_for_genome {
    my($roles,$genome,$shrub,$counts,$trip_fh) = @_;
    my @tuples = $shrub->GetAll("Genome Genome2Feature Feature Feature2Function Function Function2Role",
                                "Genome(id) = ?",
                                [$genome],
                                "Feature(id) Function(id) Function2Role(to-link)");
    my %pegs;
    foreach my $tuple (@tuples)
    {
	my($fid,$func,$role) = @$tuple;
	if ($roles->{$role} && (&SeedUtils::type_of($fid) eq 'peg'))
	{
	    $pegs{$fid}->{$role} = 1;
	}
    }
    my @fidPairs   = grep { $_->[1] } map { [$_,$shrub->loc_of($_)] } keys(%pegs);
    my @sorted     = sort { BasicLocation::Cmp($a->[1], $b->[1]) } @fidPairs;
    for (my $i=0; ($i < (@sorted-2)); $i++)
    {
	my $loc1 = $sorted[$i]->[1];
        my $loc2 = $sorted[$i+2]->[1];
        if ($loc1->Contig eq $loc2->Contig)
	{
	    &generate_triples(\%pegs,
			      $sorted[$i]->[0],
			      $sorted[$i+1]->[0],
			      $sorted[$i+2]->[0],
			      $counts,
			      $trip_fh);
	}
    }
}

sub generate_triples {
    my($pegs,$p1,$p2,$p3,$counts,$trip_fh) = @_;

    my @pairs;
    foreach my $peg ($p1,$p2,$p3)
    {
	my $h = $pegs->{$peg};
        my @roles = keys(%$h);
        push(@pairs,map { [$_,$peg] } @roles);
    }
    for (my $i=0; ($i < (@pairs-2)); $i++)
    {
	my @trp = sort { $a->[0] cmp $b->[0] } ($pairs[$i],$pairs[$i+1],$pairs[$i+2]);
	my $trip = join("\t",map { join("\t",@$_) } @trp);
        my $just_roles = join("\t",($trp[0]->[0],$trp[1]->[0],$trp[2]->[0]));
	$counts->{$just_roles}++;
	print $trip_fh $trip,"\n";
    }
}

1;
