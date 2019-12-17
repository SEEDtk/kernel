use strict;
use FIG_Config;
use SeedTkRun;

my ($dir) = @ARGV;
open (my $dh, $dir) || die "Could not open directory $dir: $!";
my @files = grep { $_ =~ /^bin\.\d+\.\d+_[12s]\.fastq/ } readdir $dh;
closedir $dh;
print scalar(@files) . " read files found.\n";
my %groups;
for my $file (@files) {
    my ($prefix, $type) = $file =~ /(bin\.\d+\.\d+)_(.)/;
    $groups{$prefix}{"-$type"} = $file;
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
    system($cmdPath, @parms);
}
