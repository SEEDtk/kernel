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
    $state->{by_vc}        = \%by_vc;
    $state->{seqs}         = \%seqs;
    $state->{to_func}      = \%to_func;
    $state->{func_to_pegs} = \%func_to_pegs;
    $state->{peg2loc}      = \%peg2loc;
    $state->{func_map}     = \%funcMap;
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
                        print "$worst: $funcName\n";
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

                $len_stats->{$func} = [ int($mean), sprintf("%0.3f",$stddev) ];
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
                my @roles_of_func = split(/[;\@\/]/, $func);
                foreach my $role_id (@roles_of_func)
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

sub project_subsys_to_genome
{
    my ( $shrub, $genome, $subsystem_id, $state, $parms ) = @_;
    print "Projecting to genome $genome.\n";
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
            $projection->{problematic_pegs}->{$peg} = $rc;
            push(@$rc,$func,$role);
        }
        else
        {
            $roles{$role}++;
            push( @$calls, [ $peg, $func, $role ] );
        }
    }
    my $pattern = &to_pattern( \%roles );
    if ($pattern)
    {

        # print STDERR "pattern=$pattern\n";
        my $poss = &possible_vc( $parms->{vc_patterns}, $pattern );

        if ($poss)
        {
            my ( $vc, $solid_genome )      = @$poss;
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
    if ( !$len ) { return ['no_translation', ''] }
    if ( $stddev < 10 ) { $stddev = 10 }
    my $delta = abs($len - $mean);
    my $max = $stddev * 3;
    if ( $delta > $max )
    {
        #print STDERR "$peg failed length test: len=$len  mean=$mean stddev=$stddev\n";
        return [ 'bad_length', "$len differs from $mean by $delta ($max = 3 * $stddev)" ];
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
      ? [ 'weak_similarity', $sim->psc . " > $worst" ]
      : ['no_similarities', ''];
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

sub all_proks {
    my($shrub) = @_;
    my @tuples = $shrub->GetAll( "Genome", "Genome(prokaryotic) = 1", [], "id" );
    my @genomes = map { $_->[0] } @tuples;
    return @genomes;
}

# returns a list of PEGs that occur within a window centered on a PEG.
# The PEG itself is returned in the list.
sub context_of_peg
{
    my ( $shrub, $peg, $window ) = @_;
    die 'context_of_peg Not yet implemented.'
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
        my $template_genome = $projection->{template_genome} || '';

        print $oh join( "\t", ( $subsystem_id, $g, $vc, $template_genome ) ), "\n";
        my $calls = $projection->{calls} || [];
        # print STDERR &Dumper($calls);
        foreach my $call ( sort { &SeedUtils::by_fig_id( $a->[0], $b->[0] ) }
            @$calls )
        {
            my ( $peg, $role, $func ) = @$call;
            print $oh "\t", join( "\t", ( $peg, $role, $func ) ), "\n";
        }
        print $oh "----\n";
        my $problems = $projection->{problematic_pegs};
        foreach
          my $peg ( sort { &SeedUtils::by_fig_id( $a, $b ) } keys(%$problems) )
        {
            print $oh join("\t",($peg,@{ $problems->{$peg} })), "\n";
        }
        print $oh "//\n";
        if ( $vc ne 'not-active' )
        {
            push @retVal, $g;
        }
    }
    return @retVal;
}

=head3 choose_genomes

    my @genomes = Projection::choose_genomes($shrub, $subsystem_id, \%badVars);

Search the database for genomes that are solid instances of the specified
subsystem. A genome is a solid instance if it has an active variant (not
C<0> or C<-1>) and all of its cells have only a single peg.

=over 4

=item shrub

The L<Shrub> object for accessing the database.

=item subsystem_id

ID of the subsytem whose genomes are desired.

=item badVars (optional)

Reference to a hash whose keys are variants to be ignored.

=item RETURN

Returns a list of genome IDs.

=back

=cut

sub choose_genomes {
    # Get the parameters.
    my ($shrub, $subsystem_id, $badVars) = @_;
    # Insure we have a bad-variant hash.
    $badVars //= {};
    # Declare the return variable.
    my @retVal;
    # Get the subsystem spreadsheet. For each cell, we want to know the genome ID and the number
    # of pegs. We filter out inactive variants.
    my @tuples = $shrub->GetAll("Subsystem2Row SubsystemRow Row2Cell Cell2Feature AND SubsystemRow Row2Genome",
                              'Subsystem2Row(from-link) = ? AND SubsystemRow(needs-curation) = ?',
                              [$subsystem_id,'0'], [qw(SubsystemRow(variant-code) Row2Genome(to-link)
                              Row2Cell(to-link) Cell2Feature(to-link))]);
    # This hash will count the number of occurrences of each subsystem cell in the list. A cell that
    # occurs more than once indicates that the genome is not solid.
    my %cells;
    # This hash will count the number of bad cells for each genome.
    my %genomes;
    # Loop through the tuples.
    for my $tuples (@tuples) {
        my ($vc, $genome, $cell, $peg) = @$tuples;
        # Only proceed if this is a good variant.
        if (! $badVars->{$vc}) {
            # Ensure every genome is represented in the hash.
            if (! exists $genomes{$genome}) {
                $genomes{$genome} = 0;
            }
            # Count this cell.
            if (++$cells{$cell} > 1) {
                $genomes{$genome}++;
            }
        }
    }
    # Return the genomes with no bad cells.
    @retVal = grep { ! $genomes{$_} } keys %genomes;
    # Return the result.
    return @retVal;
}


1;
