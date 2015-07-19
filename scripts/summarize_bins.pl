use strict;
use Data::Dumper;
use ScriptUtils;
use File::Slurp;
use Stats;

my $opt =
  ScriptUtils::Opts( '',
                     ScriptUtils::ih_options(),
                     [ 'c2uni|c=s','File connecting contigs to uni roles',{ required => 1}],
                     [ 'c1|u=f', 'incr for universal match', { default => 1 } ],
                     [ 'c2|e=f', 'decr for duplicate universal match', { default => 5} ],
                     [ 'expected=s', 'file of expected bin assignments']
    );
my $ih = ScriptUtils::IH( $opt->input );
my $c2uni = $opt->c2uni;
my $c1 = $opt->c1;
my $c2 = $opt->c2;
my $stats = Stats->new();
my %expectations;
my %expectedBin;
if ($opt->expected) {
    open(my $ih, "<", $opt->expected) || die "Could not open expectation file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        chomp $line;
        my ($contig, $bin) = split /\t/, $line;
        $expectations{$contig} = $bin;
        push @{$expectedBin{$bin}}, $contig;
    }
}
my %uni2contigs;
foreach $_ (File::Slurp::read_file($c2uni))
{
    chomp;
    my($contig,$role) = split(/\t/,$_);
    push(@{$uni2contigs{$role}},$contig);
}

$/ = "\n//\n";

while (defined ($_ = <$ih>))
{
    $/ = "\n";
    # Process the data for a bin.
    &process_1($_,\%uni2contigs);
    $/ = "\n//\n";
}
print "########\n" . $stats->Show();

sub process_1 {
    my($x,$uni2contigs) = @_;
    # Split the universal roles and the contig list.
    if ($x =~ /^(.*)\n\n(.*)\n?\/\/\n/s)
    {
        my $contigs = $1;
        my $univs = $2;
        my @univ = map { ($_ =~ /^(\d+)\t(\S.*\S)/) ? [$1,$2] : () } split(/\n/,$univs);
        my $num_univ = @univ;
        my @extra = grep { $_->[0] > 1 } @univ;
        my $num_extra = @extra;
        my %binned_contigs = map { ($_ =~ /^([^\n]+)/) ? ($1 => 1) : () } split(/\n\n/,$contigs);
        &display_univ(\@univ,$uni2contigs,\%binned_contigs);
        print "\n--------\n";
        my @counts = &count_ref($contigs);
        foreach my $tuple (@counts[0..9])
        {
            if ($tuple) {
                print join("\t",@$tuple),"\n";
            }
        }
        print "\n--------\n";
        if ($opt->expected) {
            analyze_bin([keys %binned_contigs], \%expectations, \%expectedBin);
        }
        my ($sz,$cv) = &sum_lengths($contigs);
        print "total size of contigs=$sz\n";
        print "number universal=$num_univ\n";
        print "number extra universals=$num_extra\n";
        print "average coverage=$cv\n";
         my $sc = ($c1 * $num_univ) - ($c2 * $num_extra);
        print "score=$sc\n";
        print "//\n\n";
        if ($sz >= 500000) {
            $stats->Add(Over500KBin => 1);
        } else {
            $stats->Add(Under500KBin => 1);
        }
        if ($sz >= 100000) {
            $stats->Add(Over100Kbin => 1);
        } else {
            $stats->Add(Under100Kbin => 1);
        }
    }
}

sub analyze_bin {
    my ($contigList, $expectationsH, $expectedBinHL) = @_;
    my %found;
    for my $contig (@$contigList) {
        my $bin = $expectationsH->{$contig};
        if ($bin) {
            $found{$bin}++;
        }
    }
    my @sorted = sort { $found{$b} <=> $found{$a} } keys %found;
    my $best = $sorted[0];
    if ($best) {
        my $extra = scalar(@$contigList) - $found{$best};
        print "best matching bin $best: $extra extra\n";
        for my $bin (@sorted) {
            print "* bin $bin\t$found{$bin}\n";
        }
    } else {
        print "No contigs found in expected bins.\n";
    }
}

sub display_univ {
    my($univs,$uni2contigs,$binned_contigs) = @_;
    foreach my $tuple (@$univs)
    {
        my($n,$role) = @$tuple;
        print join("\t",@$tuple),"\n";
        my $contigs = $uni2contigs->{$role};
        my %this_bin = map { ($_ => 1 ) } grep { $binned_contigs->{$_} } @$contigs;
#        if ($contigs && ($n > 1))
#        {
#            foreach my $contig (sort keys(%this_bin))
#            {
#                print "\t\t$contig\n";
#            }
#        }
    }
}

sub count_ref {
    my($contigs) = @_;

    my %counts;
    my %bp;
    my @entries = split(/\n\n/,$contigs);
    foreach my $x (@entries)
    {
        # Get the contig length and the genome name.
        my ($node, $genome) = split /\n/, $x;
        if ($node =~ /length_(\d+)/)
        {
            my $length = $1;
            my (undef, undef, $gname) = split /\t/, $genome;
            $bp{$gname} += $length;
            $counts{$gname}++;
        }
    }

    return sort { $b->[1] <=> $a->[1] } map { [$counts{$_},$bp{$_},$_] } keys(%counts);
}

sub sum_lengths {
    my($contigs) = @_;

    my $tot=0;
    my $cvg=0;
    my @entries = split(/\n\n/,$contigs);
    foreach my $x (@entries)
    {
        if ($x =~ /length_(\d+)_covg?_([\d\.]+)/)
        {
            $tot += $1;
            $cvg += $1*$2;
        }
    }
    return ($tot, $cvg/$tot);
}
