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

=head1 Analyze Occurrences of Universal Proteins

    analyze_contig_to_uni <contig.to.uni

This program analyzes the C<contig.to.uni> file produced by L<community_pipeline.pl>. It will display the number of times
each universal role was found in a contig, followed by a list of all the universal roles that were not found.

=head2 Parameters

The command-line options are those found in L<ScriptUtils/ih_options>,
which specify the full name of the C<contig.to.uni> file.

=cut

use strict;
use Stats;

my $stats = Stats->new();
my $line;
my %uni;
open(my $rh, "<$FIG_Config::global/uni.roles") || die "Could not open universal role file: $!";
while (defined($line = <$rh>)) {
    chomp $line;
    $uni{$line} = 1;
    $stats->Add($line, 0);
}
open(my $ih, "<$ARGV[0]") || die "Could not open input: $!";
while (! eof $ih) {
    $line = <$ih>; chomp $line;
    my ($contig, $uni) = split /\t/, $line;
    $stats->Add($uni, 1);
}
print $stats->Show();
print "\n-----\n";
for my $uni (sort keys %uni) {
    if (! $stats->Ask($uni)) {
        print "$uni\n";
    }
}
