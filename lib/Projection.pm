package Projection;

use gjo::stat;
use gjo::BlastInterface;
use Data::Dumper;

sub relevant_projection_data
{
    my ( $subsystem_id, $genomes, $shrub ) = @_;

    my $state = { subsystem => $subsystem_id };
    $state->{genomes} = $genomes;
    my @tuples = $shrub->GetAll(
        "Subsystem2Role Role",
        "Subsystem2Role(from-link) = ?",
        [$subsystem_id], "Role(id) Role(description)"
    );
    my %to_role_id = map { ( $_->[1] => $_->[0] ) } @tuples;
    $state->{roles} = \%to_role_id;

    @tuples = $shrub->GetAll(
        "Subsystem2Row SubsystemRow Row2Cell Cell2Feature Feature
                               Feature2Protein Protein AND
                               Feature Feature2Function Function AND
                               Feature Feature2Contig",
        "(Subsystem2Row(from-link) = ?) AND
                               (Feature2Function(security) = ?)",
        [ $subsystem_id, 2 ],
        "SubsystemRow(variant-code) Cell2Feature(to-link)
                               Function(description) Protein(sequence) Feature2Contig(to-link)
                               Feature2Contig(begin) Feature2Contig(dir)"
    );

    my %by_vc;
    my %seqs;
    my @pegs;
    my %to_func;
    my %func_to_pegs;
    my %peg2loc;

    foreach my $tuple (@tuples)
    {
        my ( $vc, $peg, $function, $seq, $contig, $begin, $strand ) = @$tuple;
        $seqs{$peg} = $seq;
        my $genome = &SeedUtils::genome_of($peg);
        $by_vc{$vc}->{$genome}->{$peg}++;
        $to_func{$peg}                   = $function;
        $func_to_pegs{$function}->{$peg} = 1;
        $peg2loc{$peg}                   = [ $contig, $begin, $strand ];
    }
    $state->{by_vc}        = \%by_vc;
    $state->{seqs}         = \%seqs;
    $state->{to_func}      = \%to_func;
    $state->{func_to_pegs} = \%func_to_pegs;
    $state->{peg2loc}      = \%peg2loc;
    return $state;
}

sub get_blast_cutoffs
{
    my ($state) = @_;
    print STDERR "computing blast cutoffs\n";

    my $func_to_pegs = $state->{func_to_pegs};
    my @funcs        = keys(%$func_to_pegs);

    # @funcs = ("Urease gamma subunit (EC 3.5.1.5)");
    my $seqs          = $state->{seqs};
    my $blast_cutoffs = {};
    foreach my $func ( sort @funcs )
    {
        print STDERR "processing function $func\n";
        my $pegH = $func_to_pegs->{$func};
        if ($pegH)
        {
            my @pegs = keys(%$pegH);
            print STDERR join( ",", @pegs ), "\n";
            if ( @pegs < 3 )    # require 3 pegs for stats
            {
                print STDERR "$func only has ", scalar @pegs, " pegs\n";
            }
            else

            {
                print STDERR "func=$func has enough pegs\n";
                my @seq_tuples;
                foreach my $peg (@pegs)
                {
                    my $tran = $seqs->{$peg};
                    if ($tran)
                    {
                        push( @seq_tuples, [ $peg, '', $tran ] );
                    }
                }
                my @output =
                  &gjo::BlastInterface::blast( \@seq_tuples, \@seq_tuples,
                    'blastp', { outForm => 'sim' } );

                my %best;
                if ( @output < 3 )
                {
                    print STDERR &Dumper( \@output, "$func has too few sims",
                        \@seq_tuples );
                }
                else
                {
                    foreach my $sim (@output)
                    {

                        # print STDERR &Dumper($sim);

                        my $id1 = $sim->id1;
                        my $id2 = $sim->id2;
                        my $sc  = $sim->psc;

                        # print STDERR "$id1 $id2 $sc\n";
                        if ( $id1 ne $id2 )
                        {
                            if (   ( !defined( $best{$id1} ) )
                                || ( $best{$id1} > $sc ) )
                            {
                                $best{$id1} = $sc;
                            }
                        }
                    }
                    my $worst = 0;
                    foreach my $id ( keys(%best) )
                    {

                      # print STDERR "processing $id best=$best{$id}  $worst\n";
                        if ( $best{$id} > $worst )
                        {
                            $worst = $best{$id};

                            # print STDERR "func=$func id=$id worst=$worst\n";
                        }
                        else
                        {

                            # print STDERR "skipping id=$id $best{$id}\n";
                        }
                    }
                    if ( $worst == 0 )
                    {
                        print STDERR "$func ha no constraints\n";
                    }
                    else
                    {
                        print STDERR "$worst: $func\n";
                    }
                    $blast_cutoffs->{$func} = [ $worst, \@seq_tuples ];
                }
            }
        }
    }
    return $blast_cutoffs;
}

sub length_stats_by_family
{
    my ($state) = @_;

    my $len_stats    = {};
    my $func_to_pegs = $state->{func_to_pegs};
    my @funcs        = keys(%$func_to_pegs);
    my $seqs         = $state->{seqs};
    foreach my $func (@funcs)
    {
        my @lengths;
        my $pegH = $func_to_pegs->{$func};
        if ($pegH)
        {
            my @pegs = keys(%$pegH);

            #           print STDERR join( ",", @pegs ), "\n";
            if ( @pegs >= 3 )    # require 3 pegs for stats
            {
                foreach my $peg (@pegs)
                {
                    my $tran = $seqs->{$peg};
                    if ($tran)
                    {
                        push( @lengths, length($tran) );
                    }
                    else
                    {
                        print STDERR "no translation for $peg\n";
                    }
                }
                my ( $mean, $stddev ) = &gjo::stat::mean_stddev(@lengths);
                $len_stats->{$func} = [ $mean, $stddev ];
                print STDERR "set mean=$mean stddev=$stddev for $func\n";
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

sub create_recognition_parameters
{
    my ( $state, $dataD, $shrub ) = @_;

    my $parms = {};
    if ( !-d $dataD )
    {
        mkdir( $dataD, 0777 ) || die "could not make $dataD";
    }
    $parms->{length_stats} = &length_stats_by_family($state);
    $parms->{blast_parms}  = &get_blast_cutoffs($state);
    $parms->{vc_patterns}  = &vc_requirements($state);

    &write_encoded_object( $parms, "$dataD/solid.projection.parms" );
}

sub vc_requirements
{
    my ($state) = @_;

    my $by_vc   = $state->{by_vc};
    my $roles   = $state->{roles};
    my $to_func = $state->{to_func};
    my %vc_patterns;
    foreach my $vc ( keys(%$by_vc) )
    {
        my $vcH     = $by_vc->{$vc};
        my @genomes = keys(%$vcH);
        foreach my $g (@genomes)
        {
            my %occ_of_role;
            my @pegs = keys( %{ $vcH->{$g} } );
            foreach my $peg (@pegs)
            {
                my $func          = $to_func->{$peg};
                my @roles_of_func = &SeedUtils::roles_of_function($func);
                foreach my $role_name (@roles_of_func)
                {
                    my $role_id;
                    if ( $role_name && ( $role_id = $roles->{$role_name} ) )
                    {
                        $occ_of_role{$role_id}++;
                        print STDERR "$g has $role_name\n";
                    }
                }
            }
            if ( keys(%occ_of_role) < 1 )
            {
                print STDERR "Genome $g has vc $vc, but no active roles\n";
            }
            else
            {
                my $pattern = &to_pattern( \%occ_of_role );
                $vc_patterns{$pattern} = [ $vc, $g ];
            }
        }
    }
    return \%vc_patterns;
}

sub from_pattern
{
    my ($pattern) = @_;

    return [ split( /,/, $pattern ) ];
}

sub to_pattern
{
    my ($roles_seen) = @_;
    return join( ",", sort keys(%$roles_seen) );
}

#########

sub write_encoded_object
{
    my ( $obj, $file ) = @_;

    my $json = JSON::XS->new;
    $json->pretty(1);
    open( OBJ, ">$file" ) || die "could not open $file";
    print OBJ $json->encode($obj);
    close(OBJ);
}

sub read_encoded_object
{
    my ($encoded_file) = @_;

    open( OBJ, "<$encoded_file" )
      || die "encoded_file $encoded_file could not be opened: $!";

    my $obj;
    my $json = JSON::XS->new;
    {
        local $/;
        undef $/;
        my $obj_txt = <OBJ>;
        $obj = $json->decode($obj_txt);
    }
    return $obj;
}

sub project_subsys_to_genome
{
    my ( $shrub, $genome, $subsystem_id, $state, $parms ) = @_;

    my $relevant           = $state->{relevant};
    my $relevant_to_genome = $relevant->{$genome};
    my @pegs               = keys(%$relevant_to_genome);

    my %roles;
    my $calls = [];
    foreach my $peg (@pegs)
    {
        my ( $role, $func, $loc ) = @{ $relevant_to_genome->{$peg} };
        if ( &good_peg( $shrub, $peg, $func, $loc, $parms ) )
        {
            $roles{$role}++;
            push( @$calls, [ $peg, $role, $func ] );
        }
    }
    my $pattern    = &to_pattern( \%roles );
    my $projection = {};
    if ($pattern)
    {
        # print STDERR "pattern=$pattern\n";
        my $poss = &possible_vc( $parms->{vc_patterns}, $pattern );
        if ($poss)
        {
            my ( $vc, $solid_genome ) = @$poss;
            $projection->{vc}              = $vc;
            $projection->{calls}           = $calls;
            $projection->{template_genome} = $solid_genome;
        }
    }
    return $projection;
}

sub possible_vc
{
    my ( $patterns, $pattern ) = @_;

    my @master_patterns = keys(%$patterns);
    my $sofar;
    my $best_pattern;

    foreach $master_pattern (@master_patterns)
    {
        if ( my $sz_master = &subset( $pattern, $master_pattern ) )
        {
            if ( ( !$sofar ) || ( $sofar > $sz_master ) )
            {
                $sofar        = $sz_master;
                $best_pattern = $master_pattern;
            }
        }
    }
    return $best_pattern ? $patterns->{$best_pattern} : undef;
}

# if set2 is a subset of set1, return size of set1; else undef
sub subset
{
    my ( $set1, $set2 ) = @_;

    my $i;
    my @set2 = split( /,/, $set2 );
    for ( $i = 0 ; ( $i < @set2 ) && ( index( $set1,$set2[$i] ) >= 0 ) ; $i++ )
    {
    }
    if ( $i == @set2 )
    {
        my $cnt = $set1 =~ tr/,/,/;
        return $cnt + 1;
    }
    return undef;
}

sub good_peg
{
    my ( $shrub, $peg, $func, $loc, $parms ) = @_;

    if ( !&ok_length( $shrub, $peg, $func, $parms ) )
    {
        print STDERR "$peg failed length check\n";
        return 0;
    }
    elsif ( !&ok_sims( $shrub, $peg, $func, $parms ) )
    {
        print STDERR "$peg failed sims check\n";
        return 0;
    }
    return 1;
}

sub ok_length
{
    my ( $shrub, $peg, $func, $parms ) = @_;

    my $tuple = $parms->{length_stats}->{$func};
    if ( !$tuple ) { return 0 }
    my ( $mean, $stddev ) = @$tuple;
    my $len = length( &seq_of_peg( $shrub, $peg ) );
    if ( !$len ) { return 0 }
    if ( $stddev < 10 ) { $stddev = 10 }
    if ( abs( ( $len - $mean ) / $stddev ) > 3 )
    {
        print STDERR
          "$peg failed length test: len=$len  mean=$mean stddev=$stddev\n";
        return 0;
    }    # z-score is too high or too low

    return 1;
}

sub ok_sims
{
    my ( $shrub, $peg, $func, $parms ) = @_;

    my $blast_parms = $parms->{blast_parms}->{$func};
    if ( !$blast_parms ) { return 1 }
    my ( $worst, $seq_tuples ) = @{$blast_parms};

    my $seq = &seq_of_peg( $shrub, $peg );
    my @sims = &gjo::BlastInterface::blast( [ $peg, '', $seq ],
        $seq_tuples, 'blastp', { outForm => 'sim' } );
    while ( ( $sim = shift @sims ) && ( $sim->id2 eq $peg ) ) { }

    # print STDERR &Dumper($sim,$sim->[10],$sim->psc);
    if ( defined($sim) && ( $sim->psc <= $worst ) )
    {
        return 1;
    }
    return 0;
}

# returns sequence(translation) of a PEG
sub seq_of_peg
{
    my ( $shrub, $peg ) = @_;

    my @tuples = $shrub->GetAll(
        "Feature2Protein Protein",
        "Feature2Protein(from-link) = ?",
        [$peg], "Protein(sequence)"
    );
    return ( @tuples > 0 ) ? $tuples[0]->[0] : undef;
}

# returns a list of PEGs that occur within a window centered on a PEG.
# The PEG itself is returned in the list.
sub context_of_peg
{
    my ( $shrub, $peg, $window ) = @_;

}

# Returns the roles of a PEG.  These are all roles from the ERDB,
# so if the PEG has extra roles (e.g., a domain with unknown function),
# you may get a reduced set compared with what you would get by just
# splitting the function.  You get back a list of 2-tuples: [id,description]
sub roles_of_peg
{
    my ( $shrub, $peg, $privilege ) = @_;

    my $roles_of_peg = [];
    my @tuples       = $shrub->GetAll(
        "Feature2Function Function Function2Role Role",
"(Feature2Function(from-link) = ?) and (Feature2Function(security) = ?)",
        [ $peg, $privilege ],
        "Function2Role(to-link) Role(description)"
    );
    return @tuples;
}

# This returns a list of roles implemented by the PEGs.
sub roles_in_peg_set
{
    my ( $shrub, $pegs, $privilege ) = @_;
    my %roles;

    foreach my $peg (@$pegs)
    {
        my @tuples = &roles_of_peg( $shrub, $peg, $privilege );
        foreach $_ (@tuples)
        {
            my ( $role_id, $desc ) = @$_;
            $roles{$desc} = $role_id;
        }
    }

    return %roles;
}

1;
