=head1 Download FASTQ Files for NCBI SRA Runs

    p3-download-samples.pl [options] outDir

This script will take as input a list of SRA sample IDs and download the FastQ files for the runs in a format
usable for the PATRIC binning service. The SRA toolkit must be installed (see L<InstallSRAtk>).

=head2 Parameters

The sole positional parameter is the name of the output directory. Each run's files will be put in a subdirectory of the
output directory having the same name as the sample ID.

The sample IDs are specified in the standard input.  The options in L<P3Utils/col_options> can be used to specify the
input column and L<P3Utils/ih_options> can be used to modify the standard input.

The following additional options are supported.

=over 4

=item missing

If specified, sample directories that already exist will be skipped.

=item clear

If specified, the output directory will be emptied before starting.

=item site

If specified, the name to be put into a C<site.tbl> file in the output folder. This should be
a lower case site name without spaces, such as is used in the HMP project.  If no site is
specified, the metadata will be read.

=item max

Maximum number of samples to download.  When this number is reached, the downloads will stop.  This
is best used with C<--missing> to download a few at a time from a specific list.

=item gzip

If specified, the output files will be in gzip format.

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
my $opt = P3Utils::script_opts('outDir', P3Utils::col_options(), P3Utils::ih_options(),
        ['missing', 'only download new samples'],
        ['clear', 'erase output directory before starting'],
        ['site=s', 'site name to specify'],
        ['gzip', 'output files should be compressed'],
        ['max=i', 'maximum number to download'],

        );
# Create a statistics object.
my $stats = Stats->new();
# Get the options.
my $missing = $opt->missing;
my $siteName;
if ($opt->site) {
    $siteName = $opt->site;
}
my $gzip = $opt->gzip;
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
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Get the sample IDs.
my $inSamples = P3Utils::get_col($ih, $keyCol);
$stats->Add(inputIDs => scalar @$inSamples);
push @samples, @$inSamples;
my $sampleTot = scalar @samples;
# Compute the number of samples to download.
my $max = $opt->max || $sampleTot;
print "$max of $sampleTot samples will be selected for downloading.\n";
my $count = 0;
# Create the SRAlib object.
my $sra = SRAlib->new(logH => \*STDOUT, stats => $stats, gzip => $gzip);
# Loop through the samples.
for my $sample (@samples) {
    $count++;
    my $target = "$outDir/$sample";
    if ($missing && -d $target) {
        print "Skipping $sample ($count of $sampleTot) with existing directory.\n";
        $stats->Add(sampleSkipped => 1);
    } elsif ($max > 0) {
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
                # Create the site file.
                my ($project, $site) = ('NCBI', $siteName);
                if (! $site) {
                    ($project, $site) = $sra->compute_site($sample);
                }
                my $siteTitle = join(' ', map { ucfirst $_ } split /_/, $site);
                # Here we need a site file.
                print "Creating site file for $siteTitle.\n";
                open(my $oh, ">$target/site.tbl") || die "Could not open site file in $target: $!";
                print $oh join("\t", $project, $site, $siteTitle) . "\n";
                # Count this download.
                $max--;
            }
        }
    }
}
print "All done:\n" . $stats->Show();
