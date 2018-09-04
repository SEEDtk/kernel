=head1 List Obsolete Genomes

    p3x-find-obsolete-genomes.pl [options] binDir

This script examines all of the binning jobs in a binning directory and lists the private genome IDs not found.
This list can be used for deletion from SOLR.

=head2 Parameters

The positional parameter is the name of the directory containing the binning sub-directories. The genome IDs will be
read from the C<index.tbl> file in each binning job's C<Eval> subdirectory.

Status messages will be written to the standard error output.

=cut

use strict;
use P3DataAPI;
use P3Utils;

# Get the command-line options.
my $opt = P3Utils::script_opts('binDir',
        );
# Get the input directory.
my ($binDir) = @ARGV;
if (! $binDir) {
    die "No input directory specified.";
} elsif (! -d $binDir) {
    die "Input directory $binDir not found or invalid.";
}
# Find all the binning jobs.
print STDERR "Scanning $binDir.\n";
opendir(my $dh, $binDir) || die "Could not open $binDir: $!";
my @bins = grep { -s "$binDir/$_/Eval/index.tbl" } readdir $dh;
closedir $dh;
print STDERR scalar(@bins) . " completed bin jobs found.\n";
# These will be the genomes to keep.
my %keep;
# Loop through the bins.
for my $bin (@bins) {
    open(my $ih, "<$binDir/$bin/Eval/index.tbl") || die "Could not open $bin: $!";
    my $count = 0;
    # Discard the header.
    my $line = <$ih>;
    # Read the genomes.
    while (! eof $ih) {
        $line = <$ih>;
        my ($sample, $id, $name) = split /\t/, $line;
        $keep{$id} = $name;
        $count++;
    }
    print STDERR "$count genomes found in $bin.\n";
}
# Get the private genomes from PATRIC.
print STDERR "Retrieving private genomes from PATRIC.\n";
my $p3 = P3DataAPI->new();
my $privates = P3Utils::get_data($p3, genome => [['eq', 'public', 0]], ['genome_id', 'genome_name']);
# Write the output headers.
P3Utils::print_cols(['genome_id', 'genome_name']);
# Output the private genomes that are not keepers.
for my $private (@$privates) {
    my ($id, $name) = @$private;
    if (! $keep{$id}) {
        P3Utils::print_cols($private);
    }
}