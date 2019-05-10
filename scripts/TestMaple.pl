use strict;
use FIG_Config;
use Stats;
use Bin::Blast;

$| = 1;
my $stats = Stats->new();
my $workDir = "$FIG_Config::temp";
open(my $ih, '<', 'empties.tbl') || die "Could not open input: $!";
while (! eof $ih) {
    my $sample = <$ih>;
    chomp $sample;
    print STDERR "Processing $sample.\n";
    $stats->Add(samplesIn => 1);
    my $blaster = Bin::Blast->new(undef, $workDir, "Bins_HMP/$sample/contigs.fasta",
            maxE => 1e-20, minlen => 0.5, gap => 600);
    # First, we need the list of bins and the locations where they hit the seed protein.
    my $matches = {};
    # We must search for the specified universal protein to create the initial bins.
    # Get the list of genomes.
    my $protFile = "$FIG_Config::global/seedprot.fa";
    $matches = $blaster->FindProtein($protFile);
    my $count = 0;
    # Save the matches to a work file.
    for my $contig (sort keys %$matches) {
        my $match = $matches->{$contig};
        print join("\t", $sample, $match->Contig, $match->Begin, $match->Dir, $match->Length) . "\n";
        $count++;
        $stats->Add('binsFound-lineOut' => 1);
    }
    print STDERR "$count seeds found for $sample.\n";
}