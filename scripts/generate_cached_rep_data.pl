use strict;
use Data::Dumper;
use gjoseqlib;
use SeedUtils;
use RepKmers;

my($complete_genomes,$sequences,$genome_index,$similarities,$Kfile);

my $usage = "usage: generate_chached_rep_data DataDir\n";
#

=head1 Generate Data for RepSeq Server

      perl generate_chached_rep_data.pl DataDir

# DataDir must exist and contain 
#
#   K                  a file containing just a number,
#                      the size of the kmers ("8" for now)         
#   complete.genomes   a 2-column table [PATRIC id,genome name]
#   6.1.1.20.fasta     a fasta file of PEG protein sequences for
#
#                      Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)
#
# This program adds
#
#   genome.index
#   similarities
#   tail.16mers
#
##########################

=cut

my $Data = shift @ARGV;
$complete_genomes = "$Data/complete.genomes";
$sequences        = "$Data/6.1.1.20.fasta";
$Kfile            = "$Data/K";
$genome_index     = "$Data/genome.index";
$similarities     = "$Data/similarities";
my $tails         = "$Data/tail.16mers";
(
 ( -d $Data)             &&
 ( -s $Kfile )           &&
 ( -s $complete_genomes) &&
 ( -s $sequences)
)
    || die $usage;

open(K,"<$Kfile") || die "could not open $Kfile";
my $K = <K>; chomp $K;
close(K);

# load vss data
open(VSS,"<$Data/vss") || die "could not open $Data/vss";
my %vss;
$/ = "\n//\n";
while (defined($_ = <VSS>))
{
    chomp;
    my @features = split(/\n/,$_);
    shift @features;
    foreach $_ (@features)
    {
	$vss{$_} = 1;
    }
}
close(VSS);
$/ = "\n";

my @seqs = grep { ! $vss{$_->[0]} } &gjoseqlib::read_fasta($sequences);
# Now convert the IDs from features to genomes
open(TAILS,">$tails") || die "could not open $tails";
foreach $_ (@seqs)
{
    $_->[0] =~ s/^fig\|(\d+\.\d+)\.peg\.\d+/$1/;
    my $tail = lc substr($_->[2],-16);
    print TAILS $_->[0],"\t",$tail,"\n";
}
close(TAILS);
my %lens = map { ($_->[0] => length($_->[2])) } @seqs;
@seqs    = sort { ($lens{$b} <=> $lens{$a}) or ($a cmp $b) } @seqs;

open(COMPLETE,"<$complete_genomes") || die "could not open $complete_genomes";
my %complete = map { (($_ =~ /^(\d+\.\d+)\s+(\S.*\S)/) && 
		      $lens{$1} && (! $vss{$1})) ? ($1 => $2) : () } <COMPLETE>;
close(COMPLETE);
my($i,$j);

my %hits_per_genome;
foreach my $tuple (@seqs)
{
    my($g,$comment,$seq) = @$tuple;
    $hits_per_genome{$g}++;
}

@seqs    = grep { $complete{$_->[0]} && ($hits_per_genome{$_->[0]} == 1) } @seqs;

my %g_to_index;
my %index_to_g;
open(INDEX,">$genome_index") || die "could not open $genome_index";
$i = 0;
foreach my $tuple (@seqs)
{
    my $g = $tuple->[0];
    if ($complete{$g} && ($hits_per_genome{$g} == 1))
    {
#	print INDEX "$i\t$g\n";
	print INDEX join("\t",($i,$g,$complete{$g})),"\n";
	$g_to_index{$g} = $i;
	$index_to_g{$i} = $g;
	$i++;
    }
}
close(INDEX);

my $Nseqs = @seqs;
my @counts;
for ($i=0; ($i < $Nseqs); $i++)
{
#    if (($i % 10) == 0) { $_ = @seqs; print STDERR "* $i $_\n" }
    for ($j=$i+1; ($j < $Nseqs); $j++)
    {
	my $n = &sim($seqs[$i]->[2],$seqs[$j]->[2],$K);
        if ($n > 3)
	{
	    push(@{$counts[$i]},[$j,$n]);
	    push(@{$counts[$j]},[$i,$n]);
	}
    }
#    print STDERR ".";
}

open(SIMS,">$similarities") || die "could not open $similarities";
$i=0;
while ($i < keys(%index_to_g))
{
    print SIMS $i,"\t",$complete{$index_to_g{$i}},"\n";
    my $counts = $counts[$i];
    if (defined($counts))
    {
	my @sorted = sort { ($b->[1] <=> $a->[1]) } @$counts;
	foreach $_ (@sorted)
	{
	    my($i2,$n) = @$_;
	    print SIMS $n,"\t",$i2,"\t",$complete{$index_to_g{$i2}},"\n";
	}
    }
    print SIMS "//\n";
    $i++;
}
close(SIMS);


sub sim {
    my($s1,$s2,$K) = @_;

    my %in1;
    my $common = 0;
    my $i;
    for ($i=0; ($i < (length($s1)-$K)); $i++)
    {
	$in1{lc substr($s1,$i,$K)} = 1;
    }

    for ($i=0; ($i < (length($s2)-$K)); $i++)
    {
	if ($in1{lc substr($s2,$i,$K)})
	{
	    $common++;
	}
    }
    return $common;
}
