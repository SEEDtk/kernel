use strict;
use Data::Dumper;
use Shrub;
use ScriptUtils;

$| = 1;
my $opt = ScriptUtils::Opts('', Shrub::script_options(),
        ['gap=i', 'region size in each direction', { default => 5000 }],
        ['roles=s', 'roles.in.subsystems file', { default => 'roles.in.subsystems'}],
        ['roleCouples=s', 'role couples file', { default => 'role.couples.tbl' }],
        ['min=i', 'minimum count', { default => 100 }]);
my $roles_in_ssF  = $opt->roles;
my $role_couplesF = $opt->rolecouples;
my $min_count     = $opt->min;
my $shrub         = Shrub->new_for_script($opt);
my $gap           = $opt->gap;

open(INSS,"<$roles_in_ssF") || die "could not open $roles_in_ssF";
my %in_subsys = map { ($_ =~ /^([^\t]+)\t[^\t]*\t(\S.*)$/) ? ($2 => $1) : () } <INSS>;
close(INSS);

my %poss;

my %roleCache;

open(RCF,"<$role_couplesF") || die "could not open $role_couplesF";
# Skip header line.
my $line = <RCF>;
# Loop through the couples.
while (defined($_ = <RCF>))
{
    chomp;
    my($r1,$r2,$count) = split(/\t/,$_);
    next if ($count < $min_count);
    next if (&ignore($r1) || &ignore($r2));

    if ($in_subsys{$r1} && (! $in_subsys{$r2}))
    {
        $poss{$r2}->{$r1} = $count;
    }
    elsif ($in_subsys{$r2} && (! $in_subsys{$r1}))
    {
        $poss{$r1}->{$r2} = $count;
    }
}

close(RCF);
my @tocheck;

foreach my $r (sort keys(%poss))
{
    push(@tocheck,[$r,$poss{$r},&best_count($poss{$r})]);
}

my @sorted = sort { $b->[2] <=> $a->[2] } @tocheck;

foreach my $tuple (@sorted)
{
    my($r1,$coupled,$count) = @$tuple;
    my $peg = &exemplar($r1,$coupled,\%roleCache,$shrub);
    print "$count\t$r1\t$peg\n";
    my @couples = sort { $coupled->{$b} <=> $coupled->{$a} } keys %$coupled;
    for my $couple (@couples) {
        print "\t$couple\t$coupled->{$couple}\n";
    }
}

sub exemplar {
    my($r1,$coupled,$roleCache,$shrub) = @_;

    my $r1_id = role_id($r1, $roleCache, $shrub);
    # The use of "2" restricts us to CoreSEED functions.
    my %pegs = map { $_ => 0 } $shrub->GetFlat('Function2Feature', 'Function2Feature(from-link) = ? AND Function2Feature(security) = ?',
            [$r1_id, 2], 'Function2Feature(to-link)');
    # Get the IDs of all the roles in which we are interested.
    my %roles = map { role_id($_, $roleCache, $shrub) => 1 } keys %$coupled;
    for my $peg (keys %pegs) {
        # Get the location of this peg.
        my $loc = $shrub->loc_of($peg);
        # Expand it by the gap distance. We might fall off the end of the contig, but we don't care, since there are no pegs there.
        $loc->Widen($gap);
        # Get all the genes in the region. The function ID is in list position 2 of each tuple.
        my $geneList = $shrub->genes_in_region($loc, 2);
        # Loop through the genes counting roles.
        for my $geneThing (@$geneList) {
            for my $role (split /[\@\;\/]/, $geneThing->[2]) {
                if ($roles{$role}) {
                    $pegs{$peg}++;
                }
            }

        }
    }
    # Find the best peg.
    my ($retVal, $count) = ('', 0);
    for my $peg (keys %pegs) {
        if ($pegs{$peg} >= $count) {
            $retVal = $peg;
            $count = $pegs{$peg};
        }
    }
    return $retVal;
}

sub role_id {
    my ($role, $roleCache, $shrub) = @_;
    my $retVal = $roleCache->{$role};
    if (! defined $retVal) {
        $retVal = $shrub->desc_to_function($role) // '';
        $roleCache->{$role} = $retVal;
    }
    return $retVal;
}

sub best_count {
    my($coupledH) = @_;

    my $best_sofar = 0;
    foreach my $close (keys(%$coupledH))
    {
        if ((! $best_sofar) || ($coupledH->{$close} > $best_sofar))
        {
            $best_sofar = $coupledH->{$close};
        }
    }
    return $best_sofar;

}

sub ignore {
    my($r) = @_;

    if ($r =~ /transport|permease|family/i)
    {
        return 1;
    }
    if ($r =~ /mobile element/i)
    {
        return 1;
    }
    return 0;
}
