use strict;
use FIG_Config;
use SeedUtils;
use File::Copy::Recursive;
use FastQ;

my $fq = FastQ->new($ARGV[0]);
my $count = 0;
while ($fq->next()) {
    $count++;
    print "$count processed.\n" if $count % 4000 == 0;
}