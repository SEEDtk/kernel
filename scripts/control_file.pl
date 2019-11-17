#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


use strict;
use warnings;
use FIG_Config;
use ScriptUtils;
use Stats;


=head1 Create a Control File for proteins.core

    control_file.pl [ options ] inDir

This script reads the training files in the specified input directory and generates a control file for creating
the associated prediction set.  The control file is tab-delimited, and contains a role ID followed by an input
width in each record.  The role ID is taken from the training file name (which is always I<roleID>C<.tbl>), and
the input width is the number of columns in the training file's header record less one (the one being the output
column).

=head2 Parameters

The positional parameter is the name of the input directory.  The control file is produced on the standard output.

The command-line options are the following

=over 4

=item verbose

Display progress messages on STDERR.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('inDir',
        ['verbose|debug|v', 'display progress on STDERR'],
        );
# Create a statistics object.
my $stats = Stats->new();
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Input directory $inDir missing or invalid.";
}
my $debug = $opt->verbose;
print STDERR "Analyzing $inDir.\n" if $debug;
opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
my @inFiles = grep { $_ =~ /^\w+\.tbl$/ } readdir $dh;
closedir $dh;
print STDERR scalar(@inFiles) . " training files found in $inDir.\n";
for my $inFile (@inFiles) {
    my ($role) = split /\./, $inFile;
    open(my $ih, '<', "$inDir/$inFile") || die "Could not open $inFile: $!";
    my $line = <$ih>;
    close $ih;
    my $cols = split(/\t/, $line) - 1;
    print "$role\t$cols\n";
    $stats->Add(roles => 1);
    $stats->Add(cols => $cols);
}
print STDERR "All done.\n" . $stats->Show() if $debug;