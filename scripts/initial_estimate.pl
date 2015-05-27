use strict;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use Shrub;
use File::Slurp;
use BasicLocation;
use Sim;
use CPscore::Vector;
use CPscore::Signal;
use CPscore::Signal::Avg;
use SampleContig;
use CPscore::Basis;
use CPscore::Basis::HotGroup;

# usage: initial_estimate -r reference.counts -d ReferenceDir > initial.estimate
my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    [ 'unicontigs|u=s','write contig to universal into this file',{}],
    [ 'refcounts|c=s', 'kmer reference counts', { required => 1 } ],
    [
        'minlen|l=i',
        'minimum length for blast match to count',
        { default => 500 }
    ],
    [
        'maxpsc|p=f',
        'maximum pscore for blast match to count',
        { default => 1.0e-100 }
    ],
    [ 'blast|b=s',  'blast type (p or n)',               { default => 'p' } ],
    [ 'minsim|s=f', 'minimum score for condensing', { default => 0.2 } ],
    [ 'savevecs|v=s', 'File used to save similarity scores',  {} ],
    [
        'refD|r=s',
        'Constructed Directory Reflecting Reference Genomes',
        { required => 1 }
    ],
    [
        'covgRatio|cr=s',
        'maximum acceptable coverage ratio for condensing',
        { default => 1.2 }
    ],
    [
        'univLimit|ul=n',
        'maximum number of duplicate universal proteins in a set', { default => 2 }
    ],
    [
        'minCovg|C=f',
        'minimum coverage amount for a community sample contig', { default => 0 }
    ],
    [ 'scoring=s', 'scoring method', { default => 'vector' }],
    [ 'compare=s', 'comparison method', { default => 'dot'}],
    [ 'logFile=s', 'name of log file'],
    [ 'blacklist=s', 'file containing reference genomes to ignore'],
    [ 'basis=s', 'method for computing the basis vectors', { default => 'normal' }],
    [ 'basisLimit=i', 'maximum number of reference genomes in the basis', { default => 3000 }],
    [ 'basisfile=s', 'output file for basis vector'],
);

my $blast_type = $opt->blast;
$blast_type = ( $blast_type =~ /^[pP]/ ) ? 'p' : 'n';

my $uni_contigsF = $opt->unicontigs;
my $ref_counts      = $opt->refcounts;
my $refD            = $opt->refd;
my $save_vecsF      = $opt->savevecs;
my $max_psc         = $opt->maxpsc;
my $min_len         = $opt->minlen;
my $min_sim         = $opt->minsim;
my $covg_constraint = $opt->covgratio;
my $univ_limit      = $opt->univlimit;
my $min_covg        = $opt->mincovg;
my $scoringType     = $opt->scoring;
my $compareType     = $opt->compare;
my $basisType       = $opt->basis;
my $blackList       = $opt->blacklist;
my $basisLimit      = $opt->basislimit;

my $scoring;
if ($scoringType eq 'vector') {
    $scoring = CPscore::Vector->new($compareType);
} elsif ($scoringType eq 'normvec') {
    $scoring = CPscore::Vector->new($compareType, normalize => 1);
} elsif ($scoringType eq 'signal') {
    $scoring = CPscore::Signal->new($compareType);
} elsif ($scoringType eq 'signalavg') {
    $scoring = CPscore::Signal::Avg->new($compareType);
} else {
    die "Invalid scoring type '$scoringType'.";
}

my $basis;
if ($basisType eq 'normal') {
    $basis = CPscore::Basis->new();
} elsif ($basisType eq 'hot') {
    $basis = CPscore::Basis::HotGroup->new();
} elsif ($basisType =~ /^hot(\d+)$/) {
    $basis = CPscore::Basis::HotGroup->new(topSize => $1);
} else {
    die "Invalid basis type '$basisType'.";
}

my %blackList;
if ($blackList) {
    %blackList = map { $_ =~ /^(\S+)/; $1 => 1 } File::Slurp::read_file($blackList, { chomp => 1 });
}
my %univ_roles = map { $_ => 1 } File::Slurp::read_file("$FIG_Config::global/uni.roles", { chomp => 1 });
print STDERR scalar(keys %univ_roles) . " universal roles read from file.\n";
opendir( REFD, $refD ) || die "Could not access $refD";
my @refs = sort { $a <=> $b } grep { $_ !~ /^\./ && ! $blackList{$_} } readdir(REFD);
my %refs = map { ( $_ => 1 ) } @refs;
closedir(REFD);

if ($opt->logfile) {
    open(LOG, ">", $opt->logfile) || die "Could not open log file: $!";
}

my @lines = File::Slurp::read_file($ref_counts);
my %ref_counts =
  map { ( ( $_ =~ /^(\d+)\t(\S+)/ ) && $refs{$1} ) ? ( $2 => $1 ) : () } @lines;
my %ref_names =
  map { ( ( $_ =~ /^\d+\t(\S+)\t(\S.*\S)/ ) && $refs{$1} ) ? ( $1 => $2 ) : () }
  @lines;
@lines = ();    # Free up memory.

my $univ_in_ref =
  &univ_roles_in_ref_pegs( $refD, \%univ_roles, \@refs, $blast_type );
my %contigs;

&process_blast_against_refs( \@refs, $refD, $univ_in_ref, $min_len, $max_psc,
    $blast_type, $min_covg, $scoring, \%contigs );

if ($uni_contigsF)
{
    &write_unis_in_contigs(\%contigs, $uni_contigsF);
}

&compute_ref_vecs( \@refs, \%contigs, $scoring, $basis, $opt->basisfile, \%ref_names, $basisLimit );
my @similarities =
  &similarities_between_contig_vecs( \%contigs, $save_vecsF, $scoring, $blast_type, $basis, $min_sim,
    $basisLimit );

my $final_sets =
  &cluster_contigs( \%contigs, \@similarities, $covg_constraint );
&output_final_sets( $final_sets, \%ref_names, \@refs, \%contigs );
if ($opt->logfile) {
    close LOG;
}

sub output_final_sets
{
    my ( $final_sets, $ref_names, $refs, $contigHash ) = @_;

    my @sets = map { $final_sets->{$_} } keys(%$final_sets);

    #   an ugly Schwarzian transform
    my @s1 =
      map { my ( $contigs, $univ ) = @$_; [ scalar keys(%$univ), $_ ] } @sets;
    my @s2 = sort { $b->[0] <=> $a->[0] } @s1;
    @sets = map { $_->[1] } @s2;

    foreach my $set (@sets)
    {
        my ( $contigs, $univ ) = @$set;
        foreach my $contig (@$contigs)
        {
            &display_contig( $contig, $ref_names, $refs,
                $contigHash->{$contig} );
        }
        &display_univ($univ);
        print "//\n";
    }
}

sub write_unis_in_contigs {
    my($contigHash, $uni_contigsF) = @_;

    open(UNI,">$uni_contigsF") || die "could not open $uni_contigsF";
    foreach my $contigID (sort keys(%$contigHash))
    {
        my $contig = $contigHash->{$contigID};
        my $x = $contig->roles;
        foreach my $univ (sort keys(%$x))
        {
            print UNI join("\t",($contig,$univ)),"\n";
        }
    }
    close(UNI);
}

sub display_univ
{
    my ($univ) = @_;

    foreach my $role ( sort keys(%$univ) )
    {
        print join( "\t", ( $univ->{$role}, $role ) ), "\n";
    }
}

sub display_contig
{
    my ( $contig, $ref_names, $refs, $contigO ) = @_;

    print "$contig\n";
    my @hits;
    my $ref_vec = $contigO->vector;
    for ( my $i = 0 ; ( $i < @$ref_vec ) ; $i++ )
    {
        if ( $ref_vec->[$i] > 0 )
        {
            push( @hits,
                [ $ref_vec->[$i], $refs->[$i], $ref_names->{ $refs->[$i] } ] );
        }
    }
    @hits = sort { ( $b->[0] <=> $a->[0] ) } @hits;
    foreach my $hit (@hits)
    {
        print join( "\t", @$hit ), "\n";
    }
    print "\n";
}

sub cluster_contigs
{
    my ( $contigHash, $similarities, $min_sim, $covg_constraint )
      = @_;

    my ( $sets, $contig_to_set ) = &initial_sets( $contigHash );
    my $final_sets =
      &condense_sets( $sets, $contig_to_set, $similarities, $min_sim,
        $covg_constraint );
    return $final_sets;
}

sub condense_sets
{
    my ( $sets, $contig_to_set, $similarities, $covg_constraint ) =
      @_;

    foreach my $sim (@$similarities)
    {
        my ( $sc, $contig1, $contig2 ) = @$sim;
        my $set1 = $contig_to_set->{$contig1};
        my $set2 = $contig_to_set->{$contig2};
        if (   $set1
            && $set2
            && ( $set1 != $set2 )
            && &univ_ok( $sets->{$set1}->[1], $sets->{$set2}->[1] )
            && &covg_ok( $sets->{$set1}[2], $sets->{$set2}[2],
                $covg_constraint ) )
        {
            if ($opt->logfile) {
                print LOG "Combining $set1 ($contig1) and $set2 ($contig2) with score $sc.\n";
            }
            my $contigs_to_move = $sets->{$set2}->[0];
            foreach my $contig_in_set2 (@$contigs_to_move)
            {
                $contig_to_set->{$contig_in_set2} = $set1;
            }
            my $contigs1 = $sets->{$set1}->[0];
            my $contigs2 = $sets->{$set2}->[0];
            push( @$contigs1, @$contigs2 );
            my $u = $sets->{$set2}->[1];
            foreach my $role ( keys(%$u) )
            {
                my $v = $sets->{$set2}->[1]->{$role};
                $sets->{$set1}->[1]->{$role} += $v;
            }

        # Compute the new coverage (covg1 * len1 + covg2 * len2) / (len1 + len2)
            $sets->{$set1}[2] =
              ( $sets->{$set1}[2] * $sets->{$set1}[3] +
                  $sets->{$set2}[2] * $sets->{$set2}[3] ) /
              ( $sets->{$set1}[3] + $sets->{$set2}[3] );
            $sets->{$set1}[3] += $sets->{$set2}[3];

            delete $sets->{$set2};
        }
    }
    return $sets;
}

sub covg_ok
{
    my ( $covg1, $covg2, $covg_constraint ) = @_;
    my $retVal;
    if ( $covg1 > $covg2 )
    {
        $retVal = ( $covg1 <= $covg_constraint * $covg2 );
    }
    elsif ( $covg1 < $covg2 )
    {
        $retVal = ( $covg2 <= $covg_constraint * $covg1 );
    }
    else
    {
        $retVal = 1;
    }
    return $retVal;
}

sub univ_ok
{
    my ( $univ1, $univ2 ) = @_;
    my $retVal = ( &in_common( $univ1, $univ2 ) <= $univ_limit );
    return $retVal;
}

sub in_common
{
    my ( $s1, $s2 ) = @_;

    my $in_common = 0;
    foreach $_ ( keys(%$s1) )
    {
        if ( $s2->{$_} )
        {
            $in_common++;
        }
    }
    return $in_common;
}

sub initial_sets
{
    my ( $contigHash ) = @_;

    my $sets          = {};
    my $contig_to_set = {};
    my $nxt_set       = 1;
    foreach my $contigID (sort keys %$contigHash)
    {
        my $contigO = $contigHash->{$contigID};
        my $univ = $contigO->roles;
        my ( $covg, $length ) = ($contigO->covg, $contigO->len);
        $sets->{$nxt_set} = [ [$contigID], $univ, $covg, $length ];
        $contig_to_set->{$contigID} = $nxt_set;
        $nxt_set++;
    }
    return ( $sets, $contig_to_set );
}

sub similarities_between_contig_vecs
{
    my ( $contigHash, $save_vecsF, $scoring, $blast_type, $basis, $min_sim, $basisLimit ) = @_;
    my @sims;
    my $typeLine = "$blast_type " . $scoring->type() . " " . $basis->type() . "-lim$basisLimit " . "$min_sim\n";
    if ( $save_vecsF && ( -s $save_vecsF ) )
    {
        open(my $ih, "<", $save_vecsF) || die "Cannot open saved vectors: $!";
        my $line = <$ih>;
        if ($line eq $typeLine) {
            @sims = sort { $b->[0] <=> $a->[0] }
          map { chomp; [ split( /\t/, $_ ) ] }
          <$ih>;
        }
    }
    if (! @sims)
    {
        @sims = &similarities_between_contig_vecs_1($contigHash, $scoring, $min_sim);
        if ( $save_vecsF && open( SAVE, ">$save_vecsF" ) )
        {
            print SAVE $typeLine;
            foreach my $tuple (@sims)
            {
                print SAVE join( "\t", @$tuple ), "\n";
            }
            close(SAVE);
        }
    }
    return @sims;
}

sub similarities_between_contig_vecs_1
{
    my ($contigHash, $scoring, $min_sim) = @_;

    my @sims;
    my @contigs = sort keys(%$contigHash);
    my $n       = @contigs;
    my ( $i, $j );
    for ( $i = 0 ; ( $i < @contigs ) ; $i++ )
    {
        my $co1 = $contigHash->{ $contigs[$i] };
        for ( $j = $i + 1 ; ( $j < @contigs ) ; $j++ )
        {
            my $co2 = $contigHash->{ $contigs[$j] };
            my $sim = $scoring->compare( $co1, $co2 );
            if ( $sim > $min_sim )
            {
                push( @sims, [ $sim, $contigs[$i], $contigs[$j] ] );
            }
        }
        if ( ( $i % 100 ) == 0 ) { print STDERR "$i of $n\n" }
    }
    print STDERR scalar(@sims) . " total scores found.\n";
    if ($scoring->sortNeeded()) {
        @sims = sort { $b->[0] <=> $a->[0] } @sims;
        print STDERR "Scores sorted.\n";
    }
    return @sims;
}

sub compute_ref_vecs
{
    my ( $refs, $contigHash, $scoring, $basis, $basisfile, $refNames, $basisLimit ) = @_;

    my $basisVec = $basis->compute($contigHash, $refs);
    if (scalar(@$basisVec) > $basisLimit) {
        my $basisLast = $basisLimit - 1;
        $basisVec = [ @{$basisVec}[0 .. $basisLast] ];
    }
    print_basis_vector($basisVec, $basisfile, $refNames);
    foreach my $contigID ( sort keys(%$contigHash) )
    {
        my $contig = $contigHash->{$contigID};
        $scoring->adjust_scores($contig);
        my $keep = $contig->Score($basisVec);
        if ( $keep )
        {
            $keep = $scoring->adjust_vector($contig);
        }
        if (! $keep)
        {
            delete $contigHash->{$contigID};
        }
    }
}

sub print_basis_vector {
    my ($basisVec, $basisfile, $refNames) = @_;
    if ($basisfile) {
       open(my $oh, ">$basisfile") || die "Could not open basis vector output file: $!";
       for my $genome (@$basisVec) {
           print $oh "$genome\t$refNames->{$genome}\n";
       }
    }
}
sub process_blast_against_refs
{
    my ( $refs, $refD, $univ_in_ref, $min_len, $max_psc, $blast_type, $min_covg, $scoring, $contigHash ) = @_;


    my $blast_out =
      ( $blast_type =~ /^[pP]/ ) ? 'blast.out.protein' : 'blast.out.dna';
    foreach my $r (sort @$refs)
    {
        my $dir = "$refD/$r";
        open( BLAST, "<$dir/$blast_out" ) || die "$dir/$blast_out is missing";
        while ( defined( $_ = <BLAST> ) )
        {
            chomp;
            my $sim = Sim->new(split( /\s+/, $_ ));
            my $contig_id = $sim->id2();
            my $contig = SampleContig::get_contig($contigHash, $contig_id);

            my $covg = $contig->covg;
            if ( ($covg >= $min_covg) && ( $sim->psc() <= $max_psc ) && ( abs( $sim->e2() - $sim->b2() ) >= $min_len ) )
            {

                $scoring->update_score($contig, $sim, $r);

                my $role;
                if ( $blast_type eq 'n' )
                {
                    if ( $role =
                        &in_univ( $univ_in_ref, $r, $sim->id1(), $sim->b1(), $sim->b2() ) )
                    {
                        $contig->SetRole($role);
                    }
                }
                else
                {
                    if ( $role = &univ_prot( $univ_in_ref, $r, $sim->id1() ) )
                    {
                        $contig->SetRole($role);
                    }
                }
            }
        }
        close(BLAST);
    }
}

sub univ_prot
{
    my ( $univ_in_ref, $ref, $ref_id ) = @_;

    if ( my $x = $univ_in_ref->{$ref}->{$ref_id} )
    {
        return $x->[0]->[1];
    }
    return undef;
}

sub in_univ
{
    my ( $univ_in_ref, $ref, $ref_contig, $beg, $end ) = @_;

    if ( my $x = $univ_in_ref->{$ref}->{$ref_contig} )
    {
        my $i;
        for (
            $i = 0 ;
            ( $i < @$x ) && ( !&between( $beg, $x->[$i]->[0], $end ) ) ;
            $i++
          )
        {
        }
        if ( $i < @$x )
        {
            return $x->[$i]->[1];
        }
    }
    return undef;
}

sub univ_roles_in_ref_pegs
{
    my ( $refD, $univ_roles, $refs, $blast_type ) = @_;

    my $univ_in_ref = {};
    foreach my $r (@$refs)
    {
        my $dir      = "$refD/$r";
        my $gto      = &SeedUtils::read_encoded_object("$dir/genome.gto");
        my $features = $gto->{features};
        foreach my $f (@$features)
        {
            if ( $univ_roles->{ $f->{function} } )
            {
                if ( $blast_type eq 'n' )
                {
                    my @locs =
                      map { BasicLocation->new($_) } @{ $f->{location} };
                    my $contig = $locs[0]->Contig;
                    my $midpt =
                      int( ( $locs[0]->Left + $locs[-1]->Right ) / 2 );
                    push(
                        @{ $univ_in_ref->{$r}->{$contig} },
                        [ $midpt, $f->{function} ]
                    );
                }
                else
                {
                    push( @{ $univ_in_ref->{$r}->{ $f->{id} } },
                        [ undef, $f->{function} ] );
                }
            }
        }
    }
    return $univ_in_ref;
}

# $univ_in_ref = &fids_of_univ_roles_in_refs( $refD, \%univ_roles, \@refs );
sub fids_of_univ_roles_in_refs {
    my ($refD, $univ_roles, $refs) = @_;

    my $univ_in_ref = {};
    foreach my $r (@$refs) {
        my $dir      = "$refD/$r";
        my $gto      = &SeedUtils::read_encoded_object("$dir/genome.gto");
        my $features = $gto->{features};
        foreach my $f (@$features)
        {
            if ( $univ_roles->{ $f->{function} } ) {
                $univ_in_ref->{$f->{id}} = $f->{function};
            }
        }
    }
    return $univ_in_ref;
}
