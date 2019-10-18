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


package Assembly;

    use strict;
    use warnings;
    use SeedUtils;
    use SeedTkRun;
    use File::Path;

=head1 Assemble Metagenomic Sample Reads

This object contains methods that assemble reads into contigs. The output is a FASTA file of the contigs and a file
describing the coverage vectors. These methods make extensive use of facilities that are currently only available on
certain Argonne unix servers.

=head2 Public Methods

=head3 Assemble

    Assembly::Assemble($outputDir, \@samples, $oh);

Assemble reads into contigs. This method assumes that an assembly rast environment is established. In particular, you
must have sourced C<~brettin/arast-user-env.sh>, used C<ar-login> to log in, and added C<~brettin/local/bin> to the path.

=over 4

=item outputDir

The output directory. The assembled contigs will be put in here in a file called C<cross.contigs>. The coverage data will
be in files called C<sampleX.bam.cov>, where I<X> is the sample number.

=item samples

Reference to a list of the names of the sample files.

=item oh (optional)

Open handle of an output file to receive status messages, or C<undef> indicating there should not be any status messages
written.

=back

=cut

sub Assemble {
    my ($outputDir, $samples, $oh) = @_;
    # First we must run the samples through the assembly RAST and get back a job ID.
    my $cmd = "ar-run -f " . join(" ", @$samples) . " -a kiki";
    my $jobID;
    my $rc = SeedTkRun::run_redirected('ar-run', '-f', @$samples, '-a', 'kiki', { stdout => \$jobID });
    if ($rc) { die "Improper exit code $rc from ar-run"; }
    chomp $jobID;
    Print($oh, "ARAST job ID is $jobID.\n");
    # Now we need to loop until the job is done.
    my ($assembled, $status);
    my $count = 10;
    while (! $assembled) {
        # Wait one minute.
        sleep 60;
        # Get the job status.
        $rc = SeedTkRun::run_redirected('ar-stat', '-j', $jobID, { stdout => \$status });
        if ($rc) { die "Improper exit code $rc from ar-stat."; }
        if ($status =~ /Complete/) {
            $assembled = 1;
        } elsif (! --$count) {
            Print($oh, "ARAST status is $status.");
            $count = 10;
        }
    }
    # Insure we have an assembly sub-directoy.
    File::Path::make_path("$outputDir/Assembly");
    # Extract the assembly.
    ##TODO
}


=head2 Internal Utility Methods

=head3 Print

    Print($oh, $text);

Write a line of text to the output stream. If the output stream is undefined, do nothing.

=over 4

=item oh

Open handle for the output stream, or C<undef> to suppress the output.

=item text

Text to write.

=back

=cut

sub Print {
    my ($oh, $text) = @_;
    if (defined $oh) {
        print $text;
    }
}

1;