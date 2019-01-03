=head1 List Good Genomes From Binning

    list_good_bins.pl [options] binDir

This script extracts a list of the good genomes from a binning directory.  This information is stored in the C<index.tbl>
file of the C<Eval> subdirectory in each binning project's directory.  The lines from this file relevant to the good genomes
(last column == 1) will be echoed to the standard output.

=head2 Parameters

The positional parameter is the name of the binning directory.

The following additional command-line options are supported.

=over 4

=item verbose

Show progress messages on the standard error output.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;


# Get the command-line options.
my $opt = P3Utils::script_opts('binDir',
        ['verbose|debug|v', 'display progress on STDERR'],
        );
# Create the statistics object.
my $stats = Stats->new();
# Get the debug flag.
my $debug = $opt->verbose;
# Get the input directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "$binDir is missing or invalid.";
}
# Get all the subdirectories.
opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
my @projects = grep { -s "$binDir/$_/Eval/index.tbl" } readdir $dh;
closedir $dh;
my $found = scalar @projects;
if (! $found) {
    die "No evaluated binning projects found in $binDir.";
} else {
    print STDERR "$found evaluated projects found in $binDir.\n" if $debug;
    # This flag will be set to TRUE after we print a header line.
    my $header;
    # This will count the projects processed.
    my $count = 0;
    # Loop through the projects.
    for my $project (@projects) {
        $count++;
        print STDERR "Processing $project ($count of $found).\n" if $debug;
        $stats->Add(projects => 1);
        # Open the index file.
        open(my $ih, '<', "$binDir/$project/Eval/index.tbl") || die "Could not open index.tbl for $project: $!";
        # Read the header line.  If this is our first time, we echo it out.
        my $line = <$ih>;
        if (! $header) {
            print $line;
            $header = 1;
        }
        # Loop through the data lines.
        while (! eof $ih) {
            $stats->Add(lineIn => 1);
            $line = <$ih>;
            if ($line =~ /1$/) {
                print $line;
                $stats->Add(lineOut => 1);
            }
        }
    }
    # All done. Output the statistics.
    print STDERR "All done.\n" . $stats->Show() if $debug;
}