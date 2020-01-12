use strict;
use FIG_Config;
use File::Copy::Recursive;
use P3DataAPI;
use P3Utils;
use GenomeTypeObject;
use GEO;

$| = 1;
my $p3 = P3DataAPI->new();
my ($dir) = @ARGV;
opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
my @files1 = readdir $dh;
closedir $dh;
my @bins;
for my $file (@files1) {
    print STDERR "$file ";
    if ($file =~ /^(\S+)\.gto$/) {
        print STDERR "kept.\n";
        push @bins, $1;
    } else {
        print STDERR "rejected.\n";
    }
}
open(my $oh, '>', "$dir/index2.tbl") || die "Could not open output file: $!";
P3Utils::print_cols(['Sample', 'Bin ID', 'Bin Name', 'Ref ID', 'Ref Name', 'Contigs', 'Base Pairs', 'N50', 'Coarse Consistency', 'Fine Consistency', 'Completeness', 'Contamination', 'Taxonomic Grouping', 'Good PheS', 'Good'], oh => $oh);
# Get the data on the genomes.
my $genomes = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name'], \@bins);
print STDERR scalar(@$genomes) . " reference genomes found.\n";
my %refs = map { $_->[0] => $_->[1] } @$genomes;
print STDERR "Extracting quality data.\n";
for my $bin (@bins) {
    my ($ref, $rName, $sample) = ('', '');
    if ($bin =~ /^\d+\.\d+/) {
        $ref = $bin;
        $rName = $refs{$bin} // '';
        $sample = "New";
    } else {
        $sample = "Old";
    }
    my $gtoName = "$dir/$bin.gto";
    if (! -s $gtoName) {
        print STDERR "$bin not found.\n";
    } else {
        print STDERR "Reading $bin.\n";
        my $gto = GenomeTypeObject->create_from_file($gtoName);
        my $geo = GEO->CreateFromGto($gtoName, detail => 0);
        my $bID = $gto->{id};
        my $bName = $gto->{scientific_name};
        my $contigs = scalar @{$gto->{contigs}};
        my $q = $gto->{quality};
        my $bp = $q->{genome_metrics}{totlen};
        my $n50 = $q->{genome_metrics}{N50};
        my $coarse = $q->{coarse_consistency};
        my $fine = $q->{fine_consistency};
        my $complt = $q->{completeness};
        my $contam = $q->{contamination};
        my $group = $q->{completeness_group} // '';
        my $goodSeed = $geo->good_seed // 0;
        my $isGood = $geo->is_good // 0;
        P3Utils::print_cols([$sample, $bID, $bName, $ref, $rName, $contigs, $bp, $n50, $coarse, $fine, $complt, $contam, $group, $goodSeed, $isGood], oh => $oh);
    }
}
