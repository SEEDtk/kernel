use strict;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use Shrub;

# usage: community_pipeline [-l MinLn] [-p MinPsc] [-s MinSim] -d DataDir -c Sample -g CalibrationGenomes
my $opt =
  ScriptUtils::Opts( '',
                     Shrub::script_options(), ScriptUtils::ih_options(),
                        [ 'minlen|l=i', 'minimum length for blast match to count', { default => 500 } ],
                        [ 'maxpsc|p=f', 'maximum pscore for blast match to count', { default => 1.0e-100 } ],
                        [ 'minsim|s=f', 'minimum % identity for condensing', { default => 7000 } ],
                        [ 'normalize', 'minimum % identity for condensing', { default => 7000 } ],
		        [ 'data|d=s', 'Data Directory for Communiy Pipeline', { required => 1 } ],
		        [ 'envsample|c=s','environmental Sample in Fasta', { } ],
		        [ 'samplename|n=s','environmental Sample Name', { required => 1 } ],
		        [ 'minkhits|k=i','minimum number hits to be a reference', { default => 400 } ],
		        [ 'refsf|r=s','File of potential reference genomes', { } ]
    );

my $dataD      = $opt->data;
my $sample     = $opt->envsample;
my $sample_id  = $opt->samplename;
my $refsF      = $opt->refsf;
my $min_hits   = $opt->minkhits;
my $max_psc    = $opt->maxpsc;
my $min_len    = $opt->minlen;
my $min_sim    = $opt->minsim;

my $ih = ScriptUtils::IH( $opt->input );
my $shrub = Shrub->new_for_script($opt);

if (! -d $dataD) { mkdir($dataD,0777) || die "could not make $dataD" }

if (! -s "$dataD/$sample_id")
{
    mkdir("$dataD/$sample_id") || die "could not make $dataD/$sample_id";
}
my $dataS = "$dataD/$sample_id";
if (! -s "$dataS/sample.fa")
{
    &SeedUtils::run("cp $sample $dataS/sample.fa");
}

if (! -d "$dataS/RefD")
{
    if (! -s $refsF)
    {
	die "you need to use the --refs parameter to specify calibration genomes";
    }
    if (! -s "$dataS/repk.json")
    {
	&SeedUtils::run("perl compute_close_data.pl > $dataS/repk.json < $refsF");
    }
    &SeedUtils::run("perl get_closest_representative_genomes.pl -d $dataS/repk.json < $dataS/sample.fa > $dataS/ref.counts");
    &SeedUtils::run("perl construct_reference_genomes.pl -c $dataS/sample.fa -m $min_hits -r $dataS/RefD < $dataS/ref.counts");
}
my $norm = $opt->normalize ? -n : '';
&SeedUtils::run("perl initial_estimate.pl -r $dataS/RefD -c $dataS/ref.counts -l $min_len -p $max_psc -s $min_sim $norm -v $dataS/saved.sim.vecs > $dataS/bins");
&SeedUtils::run("perl summarize_bins.pl < $dataS/bins > $dataS/bins.summary");
