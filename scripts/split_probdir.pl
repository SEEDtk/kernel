#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use File::Copy qw(copy);
use File::Path qw(make_path);
use SeedUtils;

my $usage = "split_probdir inDir trainDir testDir Fraction [seed]";

my $inDir;
my $trainDir;
my $testDir;
my $fraction = 0.2;
my $seed;
my $keep;
my $error;
my $help;

use Getopt::Long;
my $rc = GetOptions(
    'inDir=s'      => \$inDir,
    'trainDir=s'   => \$trainDir,
    'testDir=s'    => \$testDir,
    'fraction=f'   => \$fraction,
    'seed=i'       => \$seed,
    'keep=s'       => \$keep,
    'error=s'      => \$error,
    );
print STDERR qq(\nrc=$rc\n\n) if $ENV{VERBOSE};

if (!$rc || $help) {
    warn qq(\n   usage: $usage\n\n);
    exit(0);
}



unless ($inDir && $trainDir && $testDir && $fraction) {
    if ((@ARGV == 4) || (@ARGV == 5)) {
        ($inDir, $trainDir, $testDir, $fraction, $seed) = @ARGV;
    }
    else {
        die "Usage: $usage";
    }
}


my ($list, %keep);
if ($keep) {
    $list = 1;
    if (-s $keep) {
        %keep = map { chomp;
                      m/^(\d+)/ ? ($1 => 1) : ()
        } &SeedUtils::file_read($keep);
    }
    else {
        %keep = map { $_ => 1} split(/,/, $keep);
    }
}

if ($seed) { srand($seed); }


my ($in_X_fh, $in_y_fh, $in_ymap_fh, $in_rowh_fh);
if (!-d $inDir) {
    die "Input directory '$inDir' does not exist";
}
else {
    open($in_X_fh, "<", "$inDir/X") or die "Could not read-open '$inDir/X'";
    open($in_y_fh, "<", "$inDir/y") or die "Could not read-open '$inDir/y'";

    if (-s "$inDir/y.map") {
        open($in_ymap_fh, "<", "$inDir/y.map")
            or die "Could not read-open '$inDir/y.map'";
    }

    if (-s "$inDir/row.h") {
        open($in_rowh_fh, "<", "$inDir/row.h")
            or die "Could not read-open '$inDir/row.h'";
    }
}


my ($train_X_fh, $train_y_fh, $train_rowh_fh);
if (-d $trainDir) {
    die "Training output directory '$trainDir' already exists";
}
else {
    make_path($trainDir) || die "Could not create '$trainDir'";
    open($train_X_fh, ">", "$trainDir/X") or die "Could not write-open '$trainDir/X'";
    open($train_y_fh, ">", "$trainDir/y") or die "Could not write-open '$trainDir/y'";

    if (-s "$inDir/row.h") {
        warn "Write-opening '$trainDir/row.h'\n" if $ENV{VERBOSE};
        open($train_rowh_fh, ">", "$trainDir/row.h")
            or die "Could not write-open '$trainDir/row.h'";
    }
}


my ($test_X_fh, $test_y_fh, $test_rowh_fh);
if (-d $testDir) {
    die "test output directory '$testDir' already exists";
}
else {
    make_path($testDir) || die "Could not create '$testDir'";
    open($test_X_fh, ">", "$testDir/X") or die "Could not write-open '$testDir/X'";
    open($test_y_fh, ">", "$testDir/y") or die "Could not write-open '$testDir/y'";

    if (-s "$inDir/row.h") {
        warn "Write-opening '$testDir/row.h'\n" if $ENV{VERBOSE};
        open($test_rowh_fh, ">", "$testDir/row.h")
            or die "Could not write-open '$testDir/row.h'";
    }
}

my $line_num = 0;
my ($Xline, $yline, $ymap_line, $rowh_line);
while (defined($Xline = <$in_X_fh>) &&
       defined($yline = <$in_y_fh>)
    ) {

    my $keep;
    if ($list) {
        $keep = $keep{$line_num};
    }
    else {
        $keep = (rand() < $fraction);
    }

    if ($keep) {
        print $test_X_fh $Xline;
        print $test_y_fh $yline;

        #...NOTE: The below are dangerous, they assume a well-formatted directory
        #   with matching line-counts !!!
        if ($test_rowh_fh && defined($rowh_line = <$in_rowh_fh>)) {
            print $test_rowh_fh $rowh_line;
        }
    }
    else {
        print $train_X_fh $Xline;
        print $train_y_fh $yline;

        #...NOTE: The below are dangerous, they assume a well-formatted directory
        #   with matching line-counts !!!
        if ($train_rowh_fh && defined($rowh_line = <$in_rowh_fh>)) {
            print $train_rowh_fh $rowh_line;
        }
    }

    ++$line_num;
}


if (-s "$inDir/col.h") {
    copy("$inDir/col.h", "$testDir/col.h")
        or die "Could not copy '$inDir/col.h' to '$testDir/col.h'";

    copy("$inDir/col.h", "$trainDir/col.h")
        or die "Could not copy '$inDir/col.h' to '$trainDir/col.h'";
}


if (-s "$inDir/col.map") {
    copy("$inDir/col.map", "$testDir/col.map")
        or die "Could not copy '$inDir/col.map' to '$testDir/col.map'";

    copy("$inDir/col.map", "$trainDir/col.map")
        or die "Could not copy '$inDir/col.map' to '$trainDir/col.map'";
}


if (-s "$inDir/y.map") {
    copy("$inDir/y.map", "$testDir/y.map")
        or die "Could not copy '$inDir/y.map' to '$testDir/y.map'";

    copy("$inDir/y.map", "$trainDir/y.map")
        or die "Could not copy '$inDir/y.map' to '$trainDir/y.map'";
}

