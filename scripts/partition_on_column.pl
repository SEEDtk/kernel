#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use File::Copy qw(copy);
use File::Path qw(make_path);

use SeedUtils;

my $usage = "partition_on_column inDir outDir column";
(@ARGV == 3) || die "Usage: $usage";
my ($inDir, $outDir, $selected_column) = @ARGV;

my ($out_X_fh, $out_y_fh);
if (-d $outDir) {
    die "Output directory '$outDir' already exists";
}
else {
    make_path($outDir) || die "Could not create '$outDir'";
    open($out_X_fh, ">", "$outDir/X") or die "Could not write-open '$outDir/X'";
    open($out_y_fh, ">", "$outDir/y") or die "Could not write-open '$outDir/y'";
}


if (-s "$inDir/row.h") {
    copy("$inDir/row.h", "$outDir/row.h")
        or die "Could not copy '$inDir/row.h' to '$outDir/row.h'";
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



my ($in_X_fh, $in_y_fh, $in_ymap_fh, $in_rowh_fh);
if (!-d $inDir) {
    die "Input directory '$inDir' does not exist";
}
else {
    open($in_X_fh, "<", "$inDir/X") or die "Could not read-open '$inDir/X'";
}




my $line_num = 0;
my ($Xline, $yline, $ymap_line, $rowh_line);
while (defined($Xline = <$in_X_fh>)) {
    chomp $Xline;
    my @cells = split /\t/, $Xline;
    my $cell  = splice(@cells, $selected_column, 1);

    print $out_X_fh (join("\t", @cells), "\n");
#    print $out_y_fh ($line_num, "\t", $cell, "\n");  #...SVC script doesn't handle this format, yet... :-(
    print $out_y_fh ($cell, "\n");

    ++$line_num;
}


if (-s "$inDir/col.h") {
    my @col;
    map { m/^(\d+)\t(\S+)(\t(\S+))?/ ? ($col[$1] = [$2, ($4 || q()) ]) : () } &SeedUtils::file_read("$inDir/col.h");
    splice(@col, $selected_column, 1);

    my $colH_fh;
    open($colH_fh, '>', "$outDir/col.h")
        or die "Could not write-open '$outDir/col.h'";
    for (my $i=0; $i < @col; ++$i) {
        print $colH_fh (join("\t", ($i, $col[$i]->[0], $col[$i]->[1])), "\n");
    }
}


