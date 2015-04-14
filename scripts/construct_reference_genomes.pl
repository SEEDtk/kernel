use strict;
use Data::Dumper;
use gjoseqlib;
use SeedUtils;
use ScriptUtils;

use Shrub;

my $opt =
  ScriptUtils::Opts( '',
                     Shrub::script_options(), ScriptUtils::ih_options(),
                        [ 'contigs|c=s', 'file of sample contigs', { required => 1 } ],
		        ['minhits|m=i','minimum number of hits per ref',{ default => 400}],
                        [ 'refD|r=s', 
			  'Constructed Directory Reflecting Refernce Genomes', { required => 1 }]
    );
my $ih = ScriptUtils::IH( $opt->input );


my $contigF = $opt->contigs;
my $refD    = $opt->refd;
my $min_hits = $opt->minhits;

# usage: construct_reference_genomes -c ContigF -r ReferenceDir < close.ref.report

my @potential_orgs = grep { $_->[0] >= $min_hits } map { chop; [split(/\t/,$_)] } <$ih>;
&SeedUtils::run("formatdb -i $contigF -p F");
&pull_ref_contigs(\@potential_orgs,$contigF,$refD);

sub pull_ref_contigs {
    my($potential_orgs,$contigF,$refD) = @_;

    mkdir($refD,0777);
    foreach my $tuple (@$potential_orgs)
    {
	my($count,$g,$gs) = @$tuple;
	my $giD = "$refD/$g";
	next if (-s "$giD/blast.out");
	mkdir($giD,0777);
	&SeedUtils::run("perl /Users/rossoverbeek/Proj/SEEDtk/GIT/RASTtk/scripts/fetch_seed_gto.pl -s pseed -g $g > $giD/$g.gto");
        my $obj = &SeedUtils::read_encoded_object("$giD/$g.gto");
	my $contigs = $obj->{contigs};
	open(CONTIGS,">$giD/reference.contigs") || die "could not open $giD/reference.contigs";
	foreach my $h (@$contigs)
	{
	    print CONTIGS ">",$h->{id},"\n",$h->{dna},"\n";
	}
	close(CONTIGS);
	&SeedUtils::run("blastall -i $giD/reference.contigs -d $contigF -p blastn -FF -m 8 > $giD/blast.out");
    }
}


sub get_contigs_for_this_ref {
    my($g,$giD,$pseedO,$tmp_contigs) = @_;

    &SeedUtils::run("formatdb -i $tmp_contigs -p F");
    open(BLAST,"<$giD/blast.out") || die "could not open $giD/blast.out";
    my %keep;
    while (defined($_ = <BLAST>))
    {
	chomp;
	my @flds = split(/\s+/,$_);
	if (($flds[2] >= 70) && (abs($flds[7] - $flds[6]) > 400))
	{
	    $keep{$flds[1]} = 1;
	}
    }
    close(BLAST);

    return \%keep;
}
	
