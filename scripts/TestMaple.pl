use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use GPUtils;

print STDERR "Processing GenomePackages.\n";
my $gHash = GPUtils::get_all('GenomePackages');
print STDERR "Preparing files.\n";
open(my $ih, "<funny.tbl") || die "Could not open funny.tbl: $!";
open(my $oh, ">incorrect.tbl") || die "Could not open incorrect.tbl: $!";
my $line = <$ih>; chomp $line;
print $oh "$line\tref_genome\tlength\n";
# Read the genomes of interest.
print STDERR "Reading funny table.\n";
while (! eof $ih) {
    my ($id, @cols) = ScriptUtils::get_line($ih);
    my $gData = GPUtils::get_data($gHash, $id);
    my $refID = $gData->{'Ref Genome'};
    my $len = $gData->{'Base pairs'};
    print $oh join("\t", $id, @cols, $refID, $len) . "\n";
}
