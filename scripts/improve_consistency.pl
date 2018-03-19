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
use FastA;
use Stats;

=head1 Improve Consistency in a Genome Package

    improve_consistency.pl [ options ] packageDir

Remove contigs from a genome in order to improve its consistency. A C<contigs.tbl> report from L<analyze_consistency.pl> will
be used to identify bad contigs, which will be removed from the C<bins.fa> file. A new FASTA file will be created for submission
to RAST.

=head2 Parameters

The positional parameter is the name of the genome package directory.

The standard input file should be the C<contigs.tbl> report from L<analyze_consistency.pl>. This is a tab-delimited file with a header line
and five columns-- (0) contig ID, (1) contig length, (2) number of good features, (3) number of local features, and (4) number of bad features.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item suffix

The suffix to give to the output FASTA file. The default is C<1>, which means the output file will be C<bin1.fa> in the genome package directory.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('packageDir', ScriptUtils::ih_options(),
        ['suffix=s', 'output file name suffix', { default => '1' }]
        );
my $stats = Stats->new();
# Get the package directory.
my ($packageDir) = @ARGV;
if (! $packageDir) {
    die "No package directory specified.";
} elsif (! -d $packageDir) {
    die "Package directory $packageDir not found or invalid.";
} elsif (! -s "$packageDir/bin.fa") {
    die "$packageDir does not appear to be a genome package.";
}
# This hash will contain the IDs of the contigs to remove.
my %remove;
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
print "Reading contig file.\n";
# Throw away the header.
my $line = <$ih>;
# Loop through the contigs, finding bad ones.
while (! eof $ih) {
    my ($contig, $len, $good, $local, $bad) = ScriptUtils::get_line($ih);
    $stats->Add(contigFileIn => 1);
    if (! $good && ! $local) {
        $remove{$contig} = 1;
        $stats->Add(badContig => 1);
        $stats->Add(badLen => $len);
    } else {
        $stats->Add(goodContig => 1);
        $stats->Add(goodLen => $len);
    }
}
close $ih; undef $ih;
# Open the FASTA input file.
print "Opening FASTA files.\n";
my $fh = FastA->new("$packageDir/bin.fa");
my $suffix = $opt->suffix;
open(my $oh, ">$packageDir/bin$suffix.fa") || die "Could not open output FASTA file: $!";
print "Copying sequences.\n";
while ($fh->next) {
    $stats->Add(contigIn => 1);
    my $id = $fh->id;
    if (! $remove{$id}) {
        $fh->Write($oh);
        $stats->Add(contigOut => 1);
    }
}
print "All done:\n" . $stats->Show();
