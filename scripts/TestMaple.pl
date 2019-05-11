use strict;
use FIG_Config;
use File::Copy::Recursive;
use File::stat;
use Fcntl qw(:mode);

my ($dir) = @ARGV;
open(my $oh, '>', 'blockedFiles.tbl') || die "Could not open output: $!";
opendir(my $dh, $dir) || die "Could not open $dir: $!";
while (my $file = readdir $dh) {
    if (substr($file, 0, 1) ne '.') {
        my $path = "$dir/$file";
        my $stat = stat($path);
        if (! $stat->cando(S_IWUSR)) {
            print $oh "$path\n";
        } elsif (time - $stat->atime > 100000) {
            if (-d $path) {
                print "Erasing $path.";
                if (File::Copy::Recursive::pathempty($path)) {
                    rmdir $path;
                    print "  Done.\n";
                } else {
                    print "  Failed.\n";
                }
            } else {
                print "Deleting $path.\n";
                unlink $path;
            }
        }
    }
}