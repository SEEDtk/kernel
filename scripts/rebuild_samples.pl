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
use File::Copy::Recursive;

=head1 Copy Samples to a New Directory for ReProcessing

    rebuild_samples.pl [ options ] oldDir newDir

This is a simple script that copies basic files from an old sample directory to a new one for reprocessing. To do this, we need to transfer
the following files.

=over 4

=item 1

contigs.fasta

=item 2

output.contigs2reads.txt

=item 3

site.tbl

=item 4

XXXXXXXX_abundance_table.tsv

=back

=head2 Parameters

The positional parameters are the source and target sample directories.

The command-line options are those found in L<ScriptUtils/ih_options> plus the following.

=over 4

=item missing

Only copy samples that are missing in the output directory. The default is to erase the existing sample and re-copy.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('oldDir newDir',
        ['missing', 'only copy new samples'],
        );
my $stats = Stats->new();
# Validate the parameters.
my ($oldDir, $newDir) = @ARGV;
if (! $oldDir) {
    die "No input directory specified.";
} elsif (! -d $oldDir) {
    die "Input directory $oldDir not found or invalid.";
} elsif (! $newDir) {
    die "No output directory specified.";
} elsif (! -d $newDir) {
    die "Output directory $newDir not found or invalid.";
}
# Get the options.
my $missing = $opt->missing;
# Load the input samples.
opendir(my $dh, $oldDir) || die "Could not open $oldDir: $!";
my @samples = grep { -f "$oldDir/$_/contigs.fasta" } readdir $dh;
closedir $dh; undef $dh;
print scalar(@samples) . " input sample directories found.\n";
# Loop through the samples.
for my $sample (@samples) {
    print "Processing $sample.\n";
    $stats->Add(sampleChecked => 1);
    my $inDir = "$oldDir/$sample";
    my $outDir = "$newDir/$sample";
    my $copyOK;
    if (-d $outDir) {
        if (! $missing) {
            # Here we need to copy over an existing sample.
            $copyOK = 1;
            print "Emptying $outDir.\n";
            $stats->Add(dirCleared => 1);
            File::Copy::Recursive::pathempty($outDir) || die "Could not empty $outDir: $!";
        } else {
            print "Existing directory $outDir-- sample skipped.\n";
            $stats->Add(dirSkipped => 1);
        }
    } else {
        print "Creating $outDir.\n";
        File::Copy::Recursive::pathmk($outDir);
        $stats->Add(dirCreated => 1);
        $copyOK = 1;
    }
    if ($copyOK) {
        # Here we have a directory into which we want to copy the inputs. First we copy the main files.
        File::Copy::Recursive::fcopy("$inDir/contigs.fasta", $outDir) || die "Could not copy contig file: $!";
        File::Copy::Recursive::fcopy("$inDir/output.contigs2reads.txt", $outDir) || die "Could not copy coverage file: $!";
        print "Main files copied.\n";
        $stats->Add(samplesCopied => 1);
        if (-f "$inDir/site.tbl") {
            File::Copy::Recursive::fcopy("$inDir/site.tbl", $outDir) || die "Could not copy site file: $!";
            print "Site file copied.\n";
            $stats->Add(sitesCopied => 1);
        }
        my ($abundance) = glob("$inDir/*_abundance_table.tsv");
        if ($abundance) {
            print "Abundance file is $inDir/$abundance\n";
            File::Copy::Recursive::fcopy("$inDir/$abundance", $outDir) || die "Could not copy abundance file: $!";
            print "Abundance file copied.\n";
            $stats->Add(abundanceCopied => 1);
        }
    }
}
print "All done:\n" . $stats->Show();
