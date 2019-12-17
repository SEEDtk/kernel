use strict;
use FIG_Config;
use SeedTkRun;
use File::Copy::Recursive;

my ($dir) = @ARGV;
opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
my @files1 = readdir $dh;
closedir $dh;
my %groups;
for my $file (@files1) {
    print "$file ";
    if ($file =~ /^(bin\.\d+\.\d+)_([12s])\.fastq/) {
        print "kept.\n";
        $groups{$1}{"-$2"} = $file;
    } else {
        print "rejected.\n";
    }
}
chdir $dir;
# Find the command.
my $cmd = "spades.py";
my $cmdPath = SeedTkRun::executable_for($cmd);
die "Could not find $cmd." if ! $cmdPath;
# Windows hack for paths with spaces.
if ($cmdPath =~ /\s/) {
    $cmdPath = $cmd;
}
my @groups = sort keys %groups;
print scalar(@groups) . " bin groups found.\n";
for my $group (sort keys %groups) {
    print "Processing $group.\n";
    my $parmH = $groups{$group};
    my @parms = map { $_, $parmH->{$_} } keys %$parmH;
    File::Copy::Recursive::pathmk("asm.$group") || die "Could not make directory asm.$group: $!";
    push @parms, '-o', "asm.$group";
    print "Starting assembly for $group.\n";
    system($cmdPath, @parms);
}
