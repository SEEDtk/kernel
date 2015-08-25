#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


use strict;
use warnings;
use FIG_Config;
use ScriptUtils;
use Data::Dumper;
use SeedUtils;
my $have_pareach;
eval {
    require Proc::ParallelLoop;
    $have_pareach = 1;
};

=head1 Genetic Algorithm Parameter Search

    ga [ options ]

This script will search for the best parameters to use for a specified command using a genetic search
algorithm. The command must take as positional parameters a working directory name, an output log file
name, and one or more numeric parameters between 0 and 1. It must write to the standard output a number
indicating the degree of success (indicating the score for the specific combination of numeric parameters).
Searching for the optimal values of the numeric parameters, that is, the ones that produce the highest score,
is the purpose of this script.

=head2 Parameters

The command-line options are those found in L<ScriptUtils::ih_options> the following.

=over 4

=item command

The name of the command to execute, as described above. The command computes the score for a set of parameters.

=item kn

The number of numeric parameters. A specific set of parameter values is called a I<chromosome>. This value is
the number of values that make up a single chromosome.

=item pop_size

The number of chromosomes to generate at the beginning, called the I<starting population>. If a standard input is
specified (the C<--input> option or the default standard input), then it is presumed to be a tab-delimited file
with one chromosome per line, and each is used as a chromosome in the starting population. Chromosomes not taken
from the standard input are generated randomly. If the standard input is empty, the entire starting population is
generated randomly.

=item directory

The name of the directory to be used as input to the scoring command.

=item log

The name of a directory to contain the log files. Each log file is produced by the scoring command, and will have
a name indicating the generation number and the ID of the chromosome.

=item iterations

The number of population generations to iterate through before stopping. Each generation is created from the most
successful chromosomes of the previous generation.

=item workers

Number of parallel processes to use (Unix only).

=back

=head2 Standard Input and Output

The standard input contains starting chromosomes for the program. The file is tab-delimited, and each record should
consist of (0) the score, (1) the chromosome ID number, and (2) the one or more chromosome parameter values.

The standard output has the identical format, and contains the last generation of chromosomes. Thus, the
standard output can be used as input to a new run of the program to get improved results.

=cut

# usage: ga -p PopSz -n Interations
my $opt = ScriptUtils::Opts(
    '',
    [ 'command|c=s','Command to score a set of K-values"', { required => 1 }],
    [ 'kn|k=i','number K values in chromosome',{ default => 6 }],
    [ 'pop_size|p=i','population size',{ default => 50 }],
    [ 'directory|d=s','Directory of Precomputed Data',{ required => 1 }],
    [ 'log|l=s','Log Directory',{ required => 1 }],
    [ 'iterations|n=i','number of iterations',{ default => 100 }],
    [ 'workers|w=i', 'number of work processes', { default => 8 }]
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
    if ($have_pareach) {
        Proc::ParallelLoop::pareach(\@args,\&run,{ Max_Workers => $opt->workers });
    } else {
        for my $arg (@args) {
            &run($arg);
        }
    }

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

