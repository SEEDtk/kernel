use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use GPUtils;
use Math::Round;

my %rolesUsed;
my $totalGenomes = 0;
print STDERR "Processing GenomePackages.\n";
my $gHash = GPUtils::get_all('GenomePackages');
for my $genome (keys %$gHash) {
    my $gDir = $gHash->{$genome};
    if (-s "$gDir/EvalBySciKit/evaluate.out") {
        print STDERR "Checking $genome.\n";
        open(my $ih, "<$gDir/EvalBySciKit/evaluate.out") || die "Could not open evaluation for $genome: $!";
        while (! eof $ih) {
            my ($id, $predicted, $actual) = ScriptUtils::get_line($ih);
            if ($predicted != $actual) {
                $rolesUsed{$id}++;
            }
        }
        $totalGenomes++;
    }
}
my @roleSort = sort { $rolesUsed{$b} <=> $rolesUsed{$a} } keys %rolesUsed;
print STDERR scalar(@roleSort) . " significant roles found.\n";
print join("\t", "ID\tCount\tPercent") . "\n";
for my $role (@roleSort) {
    my $uses = $rolesUsed{$role};
    my $pct = Math::Round::nearest(0.01, $uses * 100 / $totalGenomes);
    print join("\t", $role, $uses, $pct) . "\n";
}

