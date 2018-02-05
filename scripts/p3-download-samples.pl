=head1 Download FASTQ Files for NCBI SRA Runs

    p3-download-samples.pl [options] outDir srs1 srs2 ... srsN

This script will take as input a list of SRA sample IDs and download the FastQ files for the runs in a format
usable for the PATRIC binning service. The SRA toolkit must be installed (see L<InstallSRAtk>).

=head2 Parameters

The first positional parameter is the name of the output directory. Each run's files will be put in a subdirectory of the
output directory having the same name as the sample ID.

The other positional parameters are the IDs of the samples to download. A parameter of C<-> indicates that the standard input
contains a list of the sample IDs. The options in L<P3Utils/col_options> can be used to specify the input column and
L<P3Utils/ih_options> can be used to modify the standard input.

The following additional options are supported.

=over 4

=item missing

If specified, sample directories that already exist will be skipped.

=item clear

If specified, the output directory will be emptied before starting.

=item site

If specified, the name to be put into a C<site.tbl> file in the output folder. This should be
a lower case site name without spaces, such as is used in the HMP project.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use File::Copy::Recursive;
use SRAlib;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir srs1 srs2 ... srsN', P3Utils::col_options(), P3Utils::ih_options(),
        ['missing', 'only download new samples'],
        ['clear', 'erase output directory before starting'],
        ['site=s', 'site name to specify']
        );
# Create a statistics object.
my $stats = Stats->new();
# Get the options.
my $missing = $opt->missing;
my ($siteName, $siteTitle);
if ($opt->site) {
    $siteName = $opt->site;
    $siteTitle = join(' ', map { ucfirst $_ } split /_/, $siteName);
}
# Get the output directory.
my ($outDir, @ids) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} else {
    print "Output directory is $outDir.\n";
    if ($opt->clear) {
        print "Erasing $outDir.\n";
        File::Copy::Recursive::pathempty($outDir) || die "Could not clear $outDir: $!";
    }
}
# Loop through the sample IDs.
my @samples;
for my $id (@ids) {
    if ($id eq '-') {
        # Here we want the standard input file.
        my $ih = P3Utils::ih($opt);
        # Read the incoming headers.
        my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
        my $inSamples = P3Utils::get_col($ih, $keyCol);
        $stats->Add(inputIDs => scalar @$inSamples);
        push @samples, @$inSamples;
    } else {
        push @samples, $id;
        $stats->Add(clineIDs => 1);
    }
}
my $sampleTot = scalar @samples;
print "$sampleTot samples selected for downloading.\n";
my $count = 0;
# Create the SRAlib object.
my $sra = SRAlib->new(logH => \*STDOUT, stats => $stats);
# Loop through the samples.
for my $sample (@samples) {
    $count++;
    my $target = "$outDir/$sample";
    if ($missing && -d $target) {
        print "Skipping $sample ($count of $sampleTot) with existing directory.\n";
        $stats->Add(sampleSkipped => 1);
    } else {
        print "Processing $sample ($count of $sampleTot).\n";
        my $runList = $sra->get_runs($sample);
        if (! $runList) {
            print "No runs found.\n";
            $stats->Add(sampleEmpty => 1);
        } else {
            my $runCount = scalar @$runList;
            print "Downloading $runCount runs from $sample.\n";
            my $okFlag = $sra->download_runs($runList, $target, $sample);
            if (! $okFlag) {
                print "Error downloading $sample.\n";
                $stats->Add(sampleError => 1);
            } else {
                $stats->Add(sampleDownloaded => 1);
                if ($siteName) {
                    # Here we need a site file.
                    print "Creating site file for $siteTitle.\n";
                    open(my $oh, ">$target/site.tbl") || die "Could not open site file in $target: $!";
                    print $oh join("\t", 'NCBI', $siteName, $siteTitle) . "\n";
                }
            }
        }
    }
}
print "All done:\n" . $stats->Show();
