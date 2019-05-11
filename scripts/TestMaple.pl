use strict;
use FIG_Config;
use File::Copy::Recursive;

open(my $ih, '<', 'deletes.tbl') || die "Could not open input file: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    print "Deleting $line.\n";
    File::Copy::Recursive::pathempty("Bins_HMP/$line") || die "Could not empty $line: $!";
    rmdir "Bins_HMP/$line" || die "Could not delete $line: $!";
}