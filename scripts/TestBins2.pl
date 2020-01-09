use strict;
use FIG_Config;
use SeedTkRun;
use File::Copy::Recursive;
use RASTlib;
use gjoseqlib;
use P3DataAPI;
use P3Utils;
use Bin::Improve;
use EvalHelper;

my $p3 = P3DataAPI->new();
my ($dir) = @ARGV;
opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
my @files1 = readdir $dh;
closedir $dh;
my @bins;
for my $file (@files1) {
    print "$file ";
    if ($file =~ /^asm\.bin\.(\d+\.\d+)/ && -s "$dir/$file/contigs.fasta") {
        print "kept.\n";
        push @bins, $1;
    } else {
        print "rejected.\n";
    }
}
# Get the data on the genomes.
my $genomes = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'taxon_id', 'species'], \@bins);
print scalar(@$genomes) . " genomes found.\n";
print "Initializing improver.\n";
my $improver = Bin::Improve->new($dir, p3 => $p3);
print "Processing genomes for annotation.\n";
for my $genome (@$genomes) {
    my ($bin, $taxID, $name) = @$genome;
    my $gtoName = "$dir/$bin.gto";
    if (-s $gtoName) {
        print "$bin already annotated.\n";
    } else {
        print "Reading $bin.\n";
        my $triples = gjoseqlib::read_fasta("$dir/asm.bin.$bin/contigs.fasta");
        print "Annotating $bin.\n";
        my $gto = RASTlib::Annotate($triples, $taxID, "$name clonal population", noIndex => 1);
        print "Checking eligibility for improvement.\n";
           my $ok = $improver->eligible($gto);
           if ($ok) {
               print "Improving genome.\n";
               $ok = $improver->Improve([$bin], $gto);
               if ($ok) {
                   print "Evaluating improved genome.\n";
                   EvalHelper::ProcessGto($gto, ref => $bin, p3 => $p3, external => 1, workDir => $dir);
               }
        }
        print "Writing final GTO.\n";
        SeedUtils::write_encoded_object($gto, $gtoName);
    }
}
