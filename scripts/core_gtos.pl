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
use Shrub;
use ScriptUtils;
use Shrub::GTO;
use File::Copy::Recursive;

=head1 Generate Core GTOs

    core_gtos.pl [ options ] outDir

This script fills a directory with CoreSEED L<GenomeTypeObject> files.  The L<Shrub> database is used, so the
GTOs produced are not necessarily up-to-date.

=head2 Parameters

The positional parameter is the name of the output directory.  The GTO files will be put in there.

The command-line options are those in L<Shrub/script_options> (to select the L<Shrub> database) plus
the following

=over 4

=item clear

Erase the output directory before starting.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outDir', Shrub::script_options(),
        ['clear', 'erase output directory if it exists']
        );
# Connect to the database.
print "Connecting to database.\n";
my $shrub = Shrub->new_for_script($opt);
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating output directory $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} elsif ($opt->clear) {
    print "Erasing output directory $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase $outDir: $!";
}
print "Retrieving genome IDs from Shrub.\n";
my @genomes = $shrub->GetFlat('Genome', 'Genome(core) = ? AND Genome(prokaryotic) = ?',
        [1, 1], 'id');
my $total = scalar @genomes;
my $count = 0;
for my $genome (@genomes) {
    $count++;
    print "Processing $genome ($count of $total).\n";
    my $gto = Shrub::GTO->new($shrub, $genome);
    $gto->destroy_to_file("$outDir/$genome.gto");
}
print "All done.\n";
