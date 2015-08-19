use strict;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use Proc::ParallelLoop;

# usage: ga -p PopSz -n Interations 
my $opt = ScriptUtils::Opts(
    '',
    [ 'pop_size|p=i','population size',{ default => 50 }],
    [ 'directory|d=s','Directory of Precomputed Data',{ required => 1 }],
    [ 'log|l=s','Log Directory',{ required => 1 }],
    [ 'iterations|n=i','number of iterations',{ default => 100 }]
);

my $dir        = $opt->directory;
my $logD       = $opt->log;
my $pop_sz     = $opt->pop_size;
my $iterations = $opt->iterations;
mkdir($logD,0777);

my $parms = 6;
my $pop = &init_population($parms,$pop_sz);
my $iter;
for ($iter=0; ($iter < $iterations); $iter++)
{
    $pop = &one_iteration($pop,$iter,$dir,$logD);
}
&display_final(&scored($pop,$iter,$dir,$logD));

sub display_final {
    my($scored) = @_;

    foreach my $tuple (@$scored)
    {
	print &Dumper($tuple);
    }
}

sub score {
    my($chromosome,$iter,$dir,$log,$n) = @_;
    my $args = join(" ",$dir,"$log/$iter.$n",@$chromosome);

    open(SC,"perl test_bin.pl $args |") || die "test_bin $args|";
    my $sc = <SC>;
    chomp $sc;
    close(SC);
    return $sc;;
}

#sub scored {
#    my($pop,$iter,$dir,$log) = @_;
#    my $n;
#    my @scored;
#    for ($n=0; ($n < @$pop); $n++)
#    {
#	my $chromosome = $pop->[$n];
#	my $sc = &score($chromosome,$iter,$dir,$log,$n);
#	push(@scored,[$sc,$chromosome]);
#    }
#    @scored = sort { $b->[0] <=> $a->[0] } @scored;
#    return \@scored;
#}

sub scored {
    my($pop,$iter,$dir,$log) = @_;
    my $n;
    my @args;
    my $tmpD = "TMP-$$";
    mkdir($tmpD,0777) || die "could not make $tmpD";
    for ($n=0; ($n < @$pop); $n++)
    {
	my $chromosome = $pop->[$n];
	push(@args,[$chromosome,$iter,$dir,$log,$n,$tmpD]);
    }
    pareach(\@args,\&run,{ Max_Workers => 10 });
    my @scored ;
    for (my $i=0; ($i < $n);$i++)
    {
	open(TMP,"<$tmpD/$i") || die "could not open $tmpD/$i";
	$_ = <TMP>;
	chomp;
	push(@scored,[$_,$args[$i]->[0]]);
	close(TMP);
    }
    system "rm -r $tmpD";
    @scored = sort { $b->[0] <=> $a->[0] } @scored;
    return \@scored;
}

sub run {
    my($tuple) = shift;
    my($chromosome,$iter,$dir,$log,$n,$tmpD) = @$tuple;
    my $sc = &score($chromosome,$iter,$dir,$log,$n);
    open(TMP,">$tmpD/$n") || die "could not open $tmpD/$n";
    print TMP "$sc\n";
    close(TMP);
}

sub one_iteration {
    my($pop,$iter,$dir,$logD) = @_;

    my $scored = &scored($pop,$iter,$dir,$logD);
    my @pop    = map { $_->[1] } @$scored;
    my $keep   = int(@pop / 2);
    my @new_pop;
    for (my $i=0; ($i < $keep); $i++)
    {
	push(@new_pop,$pop[$i]);
    }
    my @children;
    for (my $i=0; ($i < (@pop - $keep)); $i++)
    {
	my $p1 = int(rand() * $keep);
	my $p2 = int(rand() * $keep);
	push(@children,&mate($new_pop[$p1],$new_pop[$p2]));
    }
    push(@new_pop,@children);
    return \@new_pop;
}

sub mate {
    my($x,$y) = @_;

    my $cross = int( rand() * @$x);
    my @chrome;
    for (my $i=0; ($i < @$x); $i++)
    {
	my $v = ($i <= $cross) ? $x->[$i] : $y->[$i];
	if (rand() < 0.1)
	{
	    my $change = (rand() * 0.4) - 0.2;
	    $v      = $v + $change;
	    if ($v < 0)    { $v = 0 }
            if ($v >= $1)  { $v = 0.999 }
	}
	push(@chrome,$v);
    }
    return \@chrome;
}

sub init_population {
    my($num_parms,$pop_sz) = @_;
    
    my $pop = [];
    my $i;
    for ($i=0; ($i < $pop_sz); $i++)
    {
	my $chromosome = [];
	my $j;
	for ($j=0; ($j < $num_parms); $j++)
	{
	    my $v = rand();
	    push(@$chromosome,$v);
	}
	push(@$pop,$chromosome);
    }
    return($pop);
}

