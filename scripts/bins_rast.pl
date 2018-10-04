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
use Shrub;
use SamplePipeline;

=head1 Submit Bins Through RAST

    bins_rast.pl [ options ] workDir

This script submits a set of bins through RAST to create proposed genomes from them. The C<bins.json> file is
read to determine the bin contents, and then the C<sample.fasta> file is processed to create the FASTA files.
The C<ref.genomes.scores.tbl> is used to compute the names of the genomes (if present).

=head2 Parameters

The single positional parameter is the working directory used to access the bin data.

The command-line options are the following.

=over 4

=item user

User name for RAST access. If omitted, the default is taken from the RASTUSER environment variable.

=item password

Password for RAST access. If omitted, the default is taken from the RASTPASS environment variable.

=item sleep

Sleep interval in seconds while waiting for the job to complete. The default is C<60>.

=item partial

Only call RAST for bins for which a GTO does not yet exist. Use this to resume after a failure.

=back

=head2 Output

A JSON-format L<GenomeTypeObject> will be produced for each bin, with the name C<bin>I<X>C<.gto>, with
I<X> being the bin number.

A new copy of the bin file with universal role information embedded will be output to C<bins.rast.json>.

The analysis report will be updated in C<bins.report.txt> with the universal role data.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('workDir', Shrub::script_options(),
        ["user|u=s", "user name for RAST access"],
        ["password|p=s", "password for RAST access"],
        ["sleep=i", "sleep interval for status polling", { default => 60 }],
        ['partial', 'only process new bins'],
        );
# Verify the work directory.
my ($workDir) = @ARGV;
if (! $workDir) {
    die "No work directory specified.";
} elsif (! -d $workDir) {
    die "Invalid work directory $workDir specified.";
} else {
    my $binJsonFile = "$workDir/bins.json";
    my $contigFastaFile = "$workDir/sample.fasta";
    if (! -s $binJsonFile || ! -s $contigFastaFile) {
        die "$workDir does not appear to contain completed bins.";
    }
    # Connect to the database.
    my $shrub = Shrub->new_for_script($opt);
    # Run the RAST pipeline for bins.
    SamplePipeline::RastBins($shrub, $binJsonFile, $contigFastaFile, $workDir, user => $opt->user,
            password => $opt->password, partial => $opt->partial, sleep => $opt->sleep);
}
print "All done.\n";

