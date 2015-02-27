package Projection;

sub relevant_projection_data {
    my($subsystem_id,$genomes,$shrub) = @_;

    my $state = { subsystem => $subsystem_id };
    $state->{genomes} = $genomes;
    my @tuples = $shrub->GetAll(
	                        "Subsystem2Role Role",
	                        "Subsystem2Role(from-link) = ?",
                                [$subsystem_id], 
	                       "Role(id) Role(description)"
	                       );
    $state->{roles} = \@tuples;

    @tuples = $shrub->GetAll(
	        "Subsystem2Row SubsystemRow Row2Cell Cell2Feature Feature Feature2Protein Protein AND 
             Feature Feature2Function Function AND Feature Feature2Contig",
        "(Subsystem2Row(from-link) = ?) AND Feature2Function(security) = ?",
             [$subsystem_id,2],
        "SubsystemRow(variant-code) Cell2Feature(to-link) Function(description) Protein(sequence) 
             Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(dir)" );

    my %by_vc;
    my %seqs;
    my @pegs;
    my %to_func;
    my %func_to_pegs;;
    my %peg2loc;

    foreach my $tuple (@tuples)
    {
	my ( $vc, $peg, $function, $seq, $contig, $begin, $strand ) = @$tuple;
	$seqs{$peg} = $seq;
	my $genome = &SeedUtils::genome_of($peg);
	$by_vc{$vc}->{$genome}->{$peg}++;
	$to_func{$peg} = $function;
	$func_to_pegs{$function}->{$peg} = 1;
	$peg2loc{$peg} = [$contig,$begin,$strand];
    }
    $state->{by_vc} = \%by_vc;
    $state->{seqs}  = \%seqs;
    $state->{to_func} = \%to_func;
    $state->{func_to_pegs} = \%func_to_pegs;
    $state->{peg2loc} = \%peg2loc;
    $state->{length_stats} = &length_stats_by_family($state);
    return $state;
}
sub length_stats_by_family {
    my($state) = @_;

    my $len_stats = {};
    my $func_to_pegs = $state->{func_to_pegs};
    my @funcs        = keys(%$func_to_pegs);
    foreach my $func (@funcs)
    {
	my @lengths;
	my $pegH = $func_to_pegs->{$func};
	if ($pegH)
	{
	    my @pegs = keys(%$pegH);
	    if (@pegs >= 5)          # require 5 pegs for stats
	    {
		foreach my $peg (@pegs)
		{
		    my $tran = $seqs->{$peg};
		    if ($tran)
		    {
			push(@lengths,length($tran));
		    }
		    else
		    {
			print STDERR "no translation for $peg\n";
		    }
		}
		my($mean,$stddev) = &gjostat::mean_stddev(@lengths);
		$len_stats->{$func} = [$mean,$stddev];
	    }
	    else
	    {
		print STDERR "too few pegs for function $func\n";
	    }
	}
	else
	{
	    print STDERR "no pegs for $func\n";
	}
    }
    return $len_stats;
}

sub create_recognition_parameters {
    my($state,$dataD,$shrub) = @_;

    if (! -d $dataD)
    {
	mkdir($dataD,0777) || die "could not make $dataD";
    }
    open(PARMS,">$dataD/parameters") || die "could not open $parms";
    my $json = JSON::XS->new;
    $json->pretty(1);
    print PARMS $json->encode($state);
    close(PARMS);
}

sub project_solid_rows {
    my($state,$dataD,$shrub) = @_;

}

1;
