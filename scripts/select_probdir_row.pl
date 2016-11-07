#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use File::Copy qw(copy);
use File::Path qw(make_path);

use SeedUtils;

my $usage = "select_probdir_row inDir outDir row-index";
(@ARGV == 3) || die "Usage: $usage";
my ($inDir, $outDir, $selected_row) = @ARGV;

my $out_X_fh;
my $out_rowh_fh;
if (-d $outDir) {
    die "Output directory '$outDir' already exists";
}
else {
    make_path($outDir) || die "Could not create '$outDir'";
    open($out_X_fh,    ">", "$outDir/X")     or die "Could not write-open '$outDir/X'";
    open($out_rowh_fh, ">", "$outDir/row.h") or die "Could not write-open '$outDir/row.h'";
}


if (-s "$inDir/col.h") {
    copy("$inDir/col.h", "$outDir/col.h")
        or die "Could not copy '$inDir/col.h' to '$outDir/col.h'";
}


if (-s "$inDir/y.map") {
    copy("$inDir/y.map", "$outDir/y.map")
        or die "Could not copy '$inDir/y.map' to '$outDir/y.map'";
}


#...What was 'col.map' supposed to be used, for ???  :-(
# if (-s "$inDir/col.map") {
#     copy("$inDir/col.map", "$outDir/col.map")
# 	or die "Could not copy '$inDir/col.map' to '$outDir/col.map'";
# }



my ($in_X_fh, $in_rowh_fh);
if (!-d $inDir) {
    die "Input directory '$inDir' does not exist";
}
else {
    open($in_X_fh,    "<", "$inDir/X")     or die "Could not read-open '$inDir/X'";
    open($in_rowh_fh, "<", "$inDir/row.h") or die "Could not read-open '$inDir/row.h'";
}




my $line_num = 0;
my ($Xline, $yline, $ymap_line, $rowh_line);
while (defined($Xline = <$in_X_fh>) &&
       defined($rowh_line = <$in_rowh_fh>)
    ) {
    if (($rowh_line =~ m/^(\d+)/) && ($1 == $selected_row)) {
        print $out_X_fh    $Xline;
        print $out_rowh_fh $rowh_line;
        last;
    }
}
