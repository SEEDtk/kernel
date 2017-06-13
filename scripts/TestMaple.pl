use strict;
use FIG_Config;
use ScriptUtils;
use Shrub;
use Stats;
use SeedUtils;

my (@goodLines, @badLines);
my $stats = Stats->new();
my ($log, $dir) = @ARGV;
open(my $ih, "<$log") || die "Could not open input: $!";
while (! eof $ih) {
    my $line = <$ih>;
    if ($line =~ /(\d+) SSU ribosomal RNAs found in (\d+\.\d+)/) {
        my ($count, $genome) = ($1, $2);
        if (! open(my $qh, "<$dir/$genome/quality.tbl")) {
            print "Error opening $genome quality file: $!";
        } else {
            my $line = <$qh>;
            if ($count) {
                push @goodLines, $line;
            } else {
                push @badLines, $line;
            }
        }
    }
}
my %types = (good => \@goodLines, bad => \@badLines);
for my $type (keys %types) {
    print uc($type) . " Genome List\n";
    for my $line (@{$types{$type}}) {
        print $line;
    }
    print "\n";
}
