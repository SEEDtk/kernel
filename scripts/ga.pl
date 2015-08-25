use strict;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use Proc::ParallelLoop;

# usage: ga -p PopSz -n Interations 
my $opt = ScriptUtils::Opts(
    '',
    [ 'command|c=s','Command to score a set of K-values"',{ }],
    [ 'kn|k=i','number K values in chromosome',{ default => 6 }],
    [ 'pop_size|p=i','population size',{ default => 50 }],
    [ 'directory|d=s','Directory of Precomputed Data',{ required => 1 }],
    [ 'log|l=s','Log Directory',{ required => 1 }],
    [ 'iterations|n=i','number of iterations',{ default => 100 }]
);

my $cmd        = $opt->command;
my $parms      = $opt->kn;
my $dir        = $opt->directory;
my $logD       = $opt->log;
my $pop_sz     = $opt->pop_size;
my $iterations = $opt->iterations;
my $nxt = 1;
mkdir($logD,0777);

my $old = \*STDIN;
my $pop = &init_population($parms,$pop_sz,$cmd,\$nxt,$logD,$old);
my $iter;
for ($iter=1; ($iter <= $iterations); $iter++)
{
    $pop = &one_iteration($pop,$cmd,\$nxt,$iter,$dir,$logD);
}
&display_final($pop);

sub display_final {
    my($scored) = @_;

    foreach my $tuple (sort { $b->[0] <=> $a->[0] } @$scored)
    {
	my($sc,$id,$chromosome) = @$tuple;
	print join("\t",($sc,$id,@$chromosome)),"\n";
    }
}

sub score {
    my($cmd,$chromosome,$iter,$dir,$logD,$n) = @_;
    my $args = join(" ",$dir,"$logD/$iter.$n",@$chromosome);

    open(SC,"$cmd $args |") || die "test_bin $args|";
    my $sc = <SC>;
    chomp $sc;
    close(SC);
    return $sc;
}

#sub scored {
#    my($pop,$iter,$dir,$logD) = @_;
#    my $n;
#    my @scored;
#    for ($n=0; ($n < @$pop); $n++)
#    {
#	my $chromosome = $pop->[$n];
#	my $sc = &score($chromosome,$iter,$dir,$logD,$n);
#	push(@scored,[$sc,$chromosome]);
#    }
#    @scored = sort { $b->[0] <=> $a->[0] } @scored;
#    return \@scored;
#}

sub scored {
    my($cmd,$nxtP,$pop,$iter,$dir,$logD) = @_;
    my $n;
    my @args;
    my $tmpD = "TMP-$$";
    mkdir($tmpD,0777) || die "could not make $tmpD";
    for ($n=0; ($n < @$pop); $n++)
    {
	my $chromosome = $pop->[$n];
	push(@args,[$cmd,$chromosome,$iter,$dir,$logD,$n,$tmpD]);
    }
    pareach(\@args,\&run,{ Max_Workers => 1 });
    my @scored ;
    for (my $i=0; ($i < $n);$i++)
    {
	open(TMP,"<$tmpD/$i") || die "could not open $tmpD/$i";
	$_ = <TMP>;
	chomp;
	push(@scored,[$_,$$nxtP++,$args[$i]->[1]]);
	close(TMP);
    }
    system "rm -r $tmpD";
    @scored = sort { $b->[0] <=> $a->[0] } @scored;
    return \@scored;
}

sub run {
    my($tuple) = shift;
    my($cmd,$chromosome,$iter,$dir,$logD,$n,$tmpD) = @$tuple;
    my $sc = &score($cmd,$chromosome,$iter,$dir,$logD,$n);
    open(TMP,">$tmpD/$n") || die "could not open $tmpD/$n";
    print TMP "$sc\n";
    close(TMP);
}

sub one_iteration {
    my($pop,$cmd,$nxtP,$iter,$dir,$logD) = @_;

    my $keep   = int(@$pop / 2);
    my @new_pop;
    for (my $i=0; ($i < $keep); $i++)
    {
	push(@new_pop,$pop->[$i]);
    }
    my @children;
    for (my $i=0; ($i < (@$pop - $keep)); $i++)
    {
	my $p1 = int(rand() * $keep);
	my $p2 = int(rand() * $keep);
	push(@children,&mate($new_pop[$p1]->[2],$new_pop[$p2]->[2]));
    }
    my $children = &scored($cmd,$nxtP,\@children,$iter,$dir,$logD);
    push(@new_pop,@$children);
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
            $v = sprintf("%0.3f",$v);
	}
	push(@chrome,$v);
    }
    return \@chrome;
}

sub init_population {
    my($num_parms,$pop_sz,$cmd,$nxtP,$logD,$old) = @_;

    my @old = map { chomp; 
		    my($sc,$id,@parms) = split(/\t/,$_);
		    [$sc,$id,[@parms]] } <$old>;
    foreach $_ (@old)
    {
	if ($$nxtP <= $_->[1]) { $$nxtP = $_->[1] + 1 }
    }
    if (@old > $pop_sz) { $#old = $pop_sz -1 }

    my $pop = [];
    my $i = @old;
    while ($i < $pop_sz)
    {
	my $chromosome = [];
	my $j;
	for ($j=0; ($j < $num_parms); $j++)
	{
	    my $v = sprintf("%0.3f",rand());
	    push(@$chromosome,$v);
	}
	push(@$pop,$chromosome);
	$i++;
    }
    my $scored = &scored($cmd,$nxtP,$pop,0,$dir,$logD);
    my @sorted = sort { $b->[0] <=> $a->[0] } (@old,@$scored);
    return \@sorted;
}

