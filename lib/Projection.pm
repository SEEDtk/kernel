package Projection;

use strict;
use warnings;
use gjo::stat;
use gjo::BlastInterface;
use Data::Dumper;

sub relevant_projection_data
{
    my ( $subsystem_id, $genomes, $shrub ) = @_;

    my $state = { subsystem => $subsystem_id };
    $state->{genomes} = $genomes;
    my @tuples = $shrub->Subsystem2Role($subsystem_id);
    my %to_role_name = map { $_->[0] => $_->[1] } @tuples;
    $state->{roles} = \%to_role_name;
    my $qs = join(", ", ('?') x scalar(@$genomes));
    @tuples = $shrub->GetAll(
        "Subsystem2Row SubsystemRow Row2Cell Cell2Feature Feature
                               Feature2Protein Protein AND
                               Feature Feature2Function Function AND
                               Feature Feature2Contig AND
                               SubsystemRow Row2Genome",
        "(Subsystem2Row(from-link) = ?) AND
                               (Feature2Function(security) = ?) AND
                               Row2Genome(to-link) IN ($qs)",
        [ $subsystem_id, 2, @$genomes ],
        "SubsystemRow(variant-code) Cell2Feature(to-link)
                               Function(id) Function(description) Protein(sequence) Feature2Contig(to-link)
                               Feature2Contig(begin) Feature2Contig(dir)"
    );

    my %by_vc;
    my %seqs;
    my @pegs;
    my %to_func;
    my %func_to_pegs;
    my %peg2loc;
    my %funcMap;
    my %func_roles;

    foreach my $tuple (@tuples)
    {
        my ( $vc, $peg, $funID, $function, $seq, $contig, $begin, $strand ) = @$tuple;
        $seqs{$peg} = $seq;
        my $genome = &SeedUtils::genome_of($peg);
        $by_vc{$vc}->{$genome}->{$peg}++;
        $funcMap{$funID} = $function;
        $to_func{$peg}                   = $funID;
        $func_to_pegs{$funID}->{$peg} = 1;
        $peg2loc{$peg}                   = [ $contig, $begin, $strand ];
    }
    for my $funID (keys %funcMap) {
        $func_roles{$funID} = [ $shrub->GetFlat('Function2Role', 'Function2Role(from-link) = ?', [$funID], 'to-link') ];
    }
    $state->{by_vc}        = \%by_vc;
    $state->{seqs}         = \%seqs;
    $state->{to_func}      = \%to_func;
    $state->{func_to_pegs} = \%func_to_pegs;
    $state->{peg2loc}      = \%peg2loc;
    $state->{func_map}     = \%funcMap;
    $state->{func_roles}   = \%func_roles;
    return $state;
}

sub get_blast_cutoffs
{
    my ($state) = @_;
    #print STDERR "computing blast cutoffs\n";

    my $func_to_pegs = $state->{func_to_pegs};
    my $func_descriptions = $state->{func_map};
    my @funcs        = keys(%$func_to_pegs);

    # @funcs = ("Urease gamma subunit (EC 3.5.1.5)");
    my $seqs          = $state->{seqs};
    my $blast_cutoffs = {};
    foreach my $func ( sort @funcs )
    {
        my $funcName = "($func) $func_descriptions->{$func}";
        #print STDERR "processing function $funcName\n";
        my $pegH = $func_to_pegs->{$func};
        if ($pegH)
        {
            my @pegs = keys(%$pegH);
            #print STDERR join( ",", @pegs ), "\n";
            if ( @pegs < 3 )    # require 3 pegs for stats
            {
                print STDERR "$funcName only has ", scalar @pegs, " pegs\n";
            }
            else

            {
                #print STDERR "func=$funcName has enough pegs\n";
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
                    print STDERR &Dumper( \@output, "$funcName has too few sims",
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
                        print STDERR "$funcName has no similarity constraints\n";
                    }
                    else
                    {
                        print STDERR "$worst: $funcName\n";
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
    my $func_descriptions = $state->{func_map};
    foreach my $func (@funcs)
    {
        my $funcName = "($func) $func_descriptions->{$func}";
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
                #print STDERR "set mean=$mean stddev=$stddev for $funcName\n";
            }
            else
            {
                print STDERR "too few pegs for function $funcName\n";
            }
        }
        else
        {
            print STDERR "no pegs for $funcName\n";
        }
    }
    return $len_stats;
}

sub create_recognition_parameters
{
    my ( $state, $shrub ) = @_;

    my $parms = {};
    $parms->{length_stats} = &length_stats_by_family($state);
    $parms->{blast_parms}  = &get_blast_cutoffs($state);
    $parms->{vc_patterns}  = &vc_requirements($state);

    return $parms;
}

sub vc_requirements
{
    my ($state) = @_;
    my $by_vc   = $state->{by_vc};
    my $roles   = $state->{roles};
    my $to_func = $state->{to_func};
    my $froles  = $state->{func_roles};
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
                my $roles_of_func = $froles->{$func};
                foreach my $role_id (@$roles_of_func)
                {
                    my $role_name;
                    if ($role_name = $roles->{$role_id} )
                    {
                        $occ_of_role{$role_id}++;
                        #print STDERR "$g has $role_name\n";
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
    my ( $obj, $oh ) = @_;

# If the user passes in a file, we open it here. Because it is opened in a local
# variable, it will be closed automatically when we go out of scope. An open handle
# passed in, however, will not be closed.
    my $handle;
    if ( !ref $oh )
    {
        open( $handle, ">$oh" ) || die "Could not open output file $oh: $!";
    }
    else
    {
        $handle = $oh;
    }

    my $json = JSON::XS->new;
    $json->pretty(1);
    print $handle $json->encode($obj);

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
    print STDERR "Projecting to genome $genome.\n";
    my $relevant           = $state->{relevant};
    my $relevant_to_genome = $relevant->{$genome};
    my @pegs               = keys(%$relevant_to_genome);

    my $projection = {};
    my %roles;
    my $calls = [];
    foreach my $peg (@pegs)
    {
        my ( $role, $func, $loc ) = @{ $relevant_to_genome->{$peg} };
        my $rc = &bad_peg( $shrub, $peg, $func, $loc, $parms );
        if ($rc)
        {
            $projection->{problematic_peg}->{$peg} = $rc;
        }
        else
        {
            $roles{$role}++;
            push( @$calls, [ $peg, $role, $func ] );
        }
    }
    my $pattern = &to_pattern( \%roles );
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

    my @master_patterns = sort keys(%$patterns);
    my $sofar;
    my $best_pattern;

    foreach my $master_pattern (@master_patterns)
    {
        if ( my $sz_master = &subset( $pattern, $master_pattern ) )
        {
            if ( ( !$sofar ) || ( $sofar < $sz_master ) )
            {
                $sofar        = $sz_master;
                $best_pattern = $master_pattern;
            }
        }
    }
    return $best_pattern ? $patterns->{$best_pattern} : undef;
}

# if set2 is a subset of set1, return size of set2; else undef
sub subset
{
    my ( $set1, $set2 ) = @_;

    my $i;
    my @set2 = split( /,/, $set2 );
    for ( $i = 0 ; ( $i < @set2 ) && ( index( $set1, $set2[$i] ) >= 0 ) ; $i++ )
    {
    }
    if ( $i == @set2 )
    {
        my $cnt = $set2 =~ tr/,/,/;
        return $cnt + 1;
    }
    return undef;
}

sub bad_peg
{
    my ( $shrub, $peg, $func, $loc, $parms ) = @_;
    my $rc;
    if ( $rc = &bad_length( $shrub, $peg, $func, $parms ) )
    {
        return $rc;
    }
    elsif ( $rc = &bad_sims( $shrub, $peg, $func, $parms ) )
    {
        return $rc;
    }
    return undef;
}

sub bad_length
{
    my ( $shrub, $peg, $func, $parms ) = @_;

    my $tuple = $parms->{length_stats}->{$func};
    if ( !$tuple ) { return undef }
    my ( $mean, $stddev ) = @$tuple;
    my $len = length( &seq_of_peg( $shrub, $peg ) );
    if ( !$len ) { return ['no_translation'] }
    if ( $stddev < 10 ) { $stddev = 10 }
    if ( abs( ( $len - $mean ) / $stddev ) > 3 )
    {
        #print STDERR "$peg failed length test: len=$len  mean=$mean stddev=$stddev\n";
        return [ 'bad_length', $len, $mean, $stddev ];
    }    # z-score is too high or too low
    return undef;
}

sub bad_sims
{
    my ( $shrub, $peg, $func, $parms ) = @_;

    my $blast_parms = $parms->{blast_parms}->{$func};
    if ( !$blast_parms ) { return undef }
    my ( $worst, $seq_tuples ) = @{$blast_parms};

    my $seq = &seq_of_peg( $shrub, $peg );
    my @sims = &gjo::BlastInterface::blast( [ $peg, '', $seq ],
        $seq_tuples, 'blastp', { outForm => 'sim' } );
    my $sim;
    while ( ( $sim = shift @sims ) && ( $sim->id2 eq $peg ) ) { }

    # print STDERR &Dumper($sim,$sim->[10],$sim->psc);
    if ( defined($sim) && ( $sim->psc <= $worst ) )
    {
        return undef;
    }
    return
      defined($sim)
      ? [ 'weak_similarity', $sim->psc, $worst ]
      : ['no_similarities'];
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
    die 'context_of_peg Not yet implemented.'
}

# Returns the roles of a PEG.  These are all roles from the ERDBtk,
# so if the PEG has extra roles (e.g., a domain with unknown function),
# you may get a reduced set compared with what you would get by just
# splitting the function.  You get back a list of 2-tuples: [id,description]
sub roles_of_peg
{
    my ( $shrub, $peg, $privilege ) = @_;

    my $roles_of_peg = [];
    my @tuples       = $shrub->GetAll(
        "Feature2Function Function2Role Role",
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

sub compute_properties_of_solid_roles
{
    my ( $shrub, $subsystem_id, $genomes ) = @_;

    my $state =
      &Projection::relevant_projection_data( $subsystem_id, $genomes, $shrub );
    my $retVal = &Projection::create_recognition_parameters( $state, $shrub );
    return $retVal;
}

sub project_solid_roles
{
    my ( $shrub, $subsystem_id, $privilege, $genomes, $parms, $oh ) = @_;
    my @retVal;
    my @tuples = $shrub->GetAll(
    "Subsystem2Role Role2Function Function2Feature Feature2Contig",
        "(Subsystem2Role(from-link) = ?) AND (Function2Feature(security) = ?)",
        [ $subsystem_id, $privilege ],
        "Subsystem2Role(to-link) Subsystem2Role(ordinal) Function2Feature(to-link) Function2Feature(from-link)
                              Feature2Contig(to-link) Feature2Contig(begin) Feature2Contig(dir)"
    );
    my %funHash;
    my (%relevant, %sort);
    foreach my $tuple (@tuples)
    {
        my ( $role, $pos, $peg, $func, $contig, $beg, $strand ) = @$tuple;
        my $g = &SeedUtils::genome_of($peg);
        $sort{$role} = $pos;
        $relevant{$g}->{$peg} = [ $role, $func, [ $contig, $beg, $strand ] ];
    }

    my $state = { ( relevant => \%relevant ) };

    foreach my $g ( sort { $a <=> $b } @$genomes )
    {
        my $projection =
          &Projection::project_subsys_to_genome( $shrub, $g, $subsystem_id,
            $state, $parms );
        my $vc = $projection->{vc} ? $projection->{vc} : 'not-active';
        print $oh join( "\t", ( $subsystem_id, $g, $vc ) ), "\n";
        my $calls = $projection->{calls} || [];
        # print STDERR &Dumper($calls);
        foreach my $call ( sort { &SeedUtils::by_fig_id( $a->[0], $b->[0] ) }
            @$calls )
        {
            my ( $peg, $role, $func ) = @$call;
            print $oh "\t", join( "\t", ( $peg, $role, $func ) ), "\n";
        }
        print $oh "----\n";
        my $problems = $projection->{problematic};
        foreach
          my $peg ( sort { &SeedUtils::by_fig_id( $a, $b ) } keys(%$problems) )
        {
            print $oh @{ $problems->{$peg} }, "\n";
        }
        print $oh "//\n";
        if ( $vc ne 'not-active' )
        {
            push @retVal, $g;
        }
    }
    return @retVal;
}

1;
