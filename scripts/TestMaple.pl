use strict;
use P3DataAPI;
use P3Utils;
use EvalCom::Tax;
use GEO;
use Time::HiRes;

my $p3 = P3DataAPI->new();
my $evalCom = EvalCom::Tax->new("$FIG_Config::global/CheckG", logH => \*STDERR);
my $roleHashes = $evalCom->roleHashes;
# Get the genome IDs.
open(my $ih, '<', 'g1000.tbl') || die "Could not open input: $!";
my $header = <$ih>;
my $genomes = P3Utils::get_col($ih, 0);
my $gHash = GEO->CreateFromPatric($genomes, roleHashes => $roleHashes, detail => 0, stats => $evalCom->stats, logH => \*STDERR);
my $start = time;
my @results;
for my $genome (keys %$gHash) {
    my @stats = $evalCom->Check2($gHash->{$genome});
    push @results, [$genome, @stats];
}
my $duration = time - $start;
my $count = scalar @results;
print STDERR "$duration seconds for $count genomes.\n";
my $speed = $duration / $count;
print STDERR "$speed seconds/genome.\n";
P3Utils::print_cols(['genome', 'complete', 'contam', 'group', 'seedFlag']);
for my $result (@results) {
    P3Utils::print_cols($result);
}



