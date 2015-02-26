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
             Feature2Contig(to-link) Feature2Cointig(begin) Feature2Contig(dir)" );

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
	$by_vc{$vc}->{$peg}++;
	$to_func{$peg} = $function;
	$func_to_pegs{$function}->{$peg} = 1;
	$peg2loc{$peg) = [$contig,$begin,$strand];
    }
    $state->{by_vc} = \%by_vc;
    $state->{seqs}  = \%seqs;
    $state->{to_func} = \%to_func;
    $state->{func_to_pegs} = \%func_to_pegs;
    $state->{peg2loc} = \%peg2loc;

    return $state;
}

sub create_recognition_parameters {
    my($state,$dataD,$shrub) = @_;

}

sub project_solid_rows {
    my($state,$dataD,$shrub) = @_;

}

1;
