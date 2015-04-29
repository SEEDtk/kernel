use strict;
use Data::Dumper;
use ScriptUtils;

my $opt =
  ScriptUtils::Opts( '',
                     ScriptUtils::ih_options(),
                     [ 'c1|u=f', 'incr for universal match', { default => 1 } ],
                     [ 'c2|e=f', 'decr for duplicate universal match', { default => 5} ]
    );
my $ih = ScriptUtils::IH( $opt->input );
my $c1 = $opt->c1;
my $c2 = $opt->c2;

$/ = "\n//\n";

while (defined ($_ = <STDIN>))
{
    $/ = "\n";
    # Process the data for a bin.
    &process_1($_);
    $/ = "\n//\n";
}

sub process_1 {
    my($x) = @_;
    # Split the universal roles and the contig list.
    if ($x =~ /^(.*)\n\n(.*)\n?\/\/\n/s)
    {
        my $contigs = $1;
        my $univs = $2;
	my @univ = map { ($_ =~ /^(\d+)/) ? $1 : () } split(/\n/,$univs);
	my $num_univ = @univ;
	my @extra = grep { $_ > 1 } @univ;
	my $num_extra = @extra;
        print $univs,"\n--------\n";
        my @counts = &count_ref($contigs);
        foreach my $tuple (@counts)
        {
            print join("\t",@$tuple),"\n";
        }
        print "\n--------\n";
        my $sz = &sum_lengths($contigs);
        print "total size of contigs=$sz\n";
	print "number universal=$num_univ\n";
	print "number extra universals=$num_extra\n";
	my $sc = ($c1 * $num_univ) - ($c2 * $num_extra);
	print "score=$sc\n";
        print "//\n\n";
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
        if ($x =~ /length_(\d+).*\t([^\t]*)$/s)
        {
            $bp{$2} += $1;
            $counts{$2}++;
        }
    }

    return sort { $b->[0] <=> $a->[0] } map { [$counts{$_},$bp{$_},$_] } keys(%counts);
}

sub sum_lengths {
    my($contigs) = @_;

    my $tot=0;
    my @entries = split(/\n\n/,$contigs);
    foreach my $x (@entries)
    {
        if ($x =~ /length_(\d+)/)
        {
            $tot += $1;
        }
    }
    return $tot;
}
