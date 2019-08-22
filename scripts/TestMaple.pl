use strict;
use FIG_Config;
use SeedUtils;
use P3Utils;
use GEO;
use EvalCon;
use Stats;

my $stats = Stats->new();
print STDERR "Loading role hashes.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
print STDERR "Loading GTO directory.\n";
opendir(my $dh, "Italians/0_10000") || die "Could not open italian directory: $!";
open(my $ih, "<italians2.tbl") || die "Could not open input: $!";
my @gtos = grep { $_ =~ /^\d+\.\d+\.gto$/ } readdir $dh;
my %gtos;
for my $gto (@gtos) {
    if ($gto =~ /(\d+\.\d+)\.gto/) {
        $gtos{$1} = $gto;
    }
}
my %geoOptions = (roleHashes => [$nMap, $cMap], stats => $stats, detail => 0, logH => \*STDERR);
# print the output header
my @roles = sort keys %$nMap;
P3Utils::p3_cols(['genome_id', 'genome_name', 'completeness', 'contamination', @roles]);
print STDERR scalar(@roles) . " roles found in database.\n";
my ($headers, $cols) = P3Utils::find_headers($ih, input => 'genome_id', 'genome_name', 'completeness', 'contamination');
while (! eof $ih) {
    my ($genome, $name, $complete, $contam) = P3Utils::get_cols($ih, $cols);
    $stats->Add(genomeIn => 1);
    if ($gtos{$genome}) {
        print STDERR "Processing $genome: $name\n";
        $stats->Add(gtoFound => 1);
        my $gto = GenomeTypeObject->create_from_file($gtos{$genome});
        my $geo = GEO->CreateFromGto($gto, %geoOptions);
        my $roleH = $geo->roleCounts;
        P3Utils::p3_cols([$genome, $name, $complete, $contam, map { $roleH->{$_} // 0 } @roles]);
    }
}
print STDERR "All done.\n" . $stats->Show();
