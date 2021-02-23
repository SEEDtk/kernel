use strict;
use FIG_Config;
use File::Copy::Recursive;

# parms are dir1 dir2 dirO

# Tracks files already in the output
my $outH;
my ($dir1, $dir2, $dirO) = @ARGV;
# Verify the output directory.
if (! -d $dirO) {
    File::Copy::Recursive::pathmk($dirO) || die "Could not create $dirO: $!";
} else {
    $outH = getFastQ($dirO);
}
# Get the two input directories.
my $dir1H = getFastQ($dir1);
my $dir2H = getFastQ($dir2);
for my $file1 (sort keys %$dir1H) {
    my $file2 = $file1;
    $file2 =~ s/_L001/_L002/;
    if ($dir2H->{$file2}) {
        my $file3 = $file1;
        $file3 =~ s/_L001/_L003/;
        if (-s "$dirO/$file3") {
            print "Skipping $file3-- already done.\n";
        } else {
            print "Copying $file1 and $file2 to $file3.\n";
            open(my $oh, '>', "$dirO/$file3") || die "Could not open output file $file3: $!";
            copyFile("$dir1/$file1", $oh);
            copyFile("$dir2/$file2", $oh);
            close $oh;
        }
    }
}
print "All done.\n";

# Copy a file to the output.
sub copyFile {
    my ($file, $oh) = @_;
    open(my $ih, '<', $file) || die "Could not open input file $file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        print $oh $line;
    }
}


# Return all the FASTQ files in the directory as a hash.
sub getFastQ {
    my ($dir) = @_;
    opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
    my @files = grep { $_ =~ /\.fastq$/ } readdir $dh;
    my %retVal = map { $_ => 1 } @files;
    return \%retVal;
}
