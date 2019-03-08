=head1 Find Good Genomes After Downloadable-Genome Processing

    p3x-find-processed-genomes.pl [options] mainDir

This script takes the output file from L<p3x-process-downloadable-genomes.pl> and creates a list of all the good genomes along
with an indication of where to find the L<GenomeTypeObject> and FASTA files produced during the run.  This list can then be
used to load the genomes into PATRIC.

=head2 Parameters

The single positional parameter is the name of the master output directory.  The GTO and FASTA files should all be in this
directory or one of its subdirectories.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The additional command-line options are the following.

=over 4

=item verbose

Write progress messages to STDERR.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use GEO;

# Get the command-line options.
my $opt = P3Utils::script_opts('mainDir', P3Utils::ih_options(),
        ['verbose|debug|v', 'write progress messages to STDERR']
        );
# Get the debug option.
my $debug = $opt->verbose;
my $stats = Stats->new();
# Get the main directory.
my ($mainDir) = @ARGV;
if (! $mainDir) {
    die "No main directory specified.";
} elsif (! -d $mainDir) {
    die "Invalid or missing directory $mainDir: $!";
}
# Get a directory of all the FASTA and GTO files.
print STDERR "Analyzing input directory $mainDir.\n" if $debug;
my (%fastaFiles, %gtoFiles);
my @dirStack = ($mainDir);
while (@dirStack) {
    my $dir = pop @dirStack;
    $stats->Add(dirIn => 1);
    opendir(my $dh, $dir) || die "Could not open $dir: $!";
    while (my $file = readdir $dh) {
        if (substr($file, 0, 1) ne '.') {
            $stats->Add(fileIn => 1);
            # Form the full file name.
            my $fileName = "$dir/$file";
            if (-d $fileName) {
                # Another directory, so stack it.
                push @dirStack, $fileName;
            } elsif ($file =~ /^(.+)\.fa$/) {
                # FASTA file, goes in the FASTA hash.
                $fastaFiles{$1} = $fileName;
                $stats->Add(fastaIn => 1);
            } elsif ($file =~ /^(.+)\.gto$/) {
                # GTO file, goes in the GTO hash.
                $gtoFiles{$1} = $fileName;
                $stats->Add(gtoIn => 1);
            }
        }
    }
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# We need a hash of all the genome IDs for genomes already in PATRIC.
my $gHash = get_all_genomes($p3);
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $cols) = P3Utils::find_headers($ih, input => 'genome_name', 'PATRIC ID', 'Scientific Name', 'Good?');
# Write output headers.
P3Utils::print_cols(['genome_id', 'genome_name', 'FASTA', 'GTO', 'PATRIC_status']);
# Loop through the input.
print STDERR "Processing input.\n" if $debug;
my $count = 0;
while (! eof $ih) {
    # Get the next genome.
    my ($gName, $genomeID, $genomeName, $goodFlag) = P3Utils::get_cols($ih, $cols);
    $stats->Add(genomeIn => 1);
    # Only proceed if it's good.
    if (! $goodFlag) {
        $stats->Add(genomeBadScore => 1);
    } else {
        # Compute the PATRIC status.
        my $pStatus = $gHash->{$genomeID};
        if ($pStatus) {
            $stats->Add(genomePatric => 1);
        } else {
            $pStatus = '';
        }
        # Get the GTO.
        my $gtoFile = $gtoFiles{$genomeID} // '';
        if (! $gtoFile) {
            $stats->Add(genomeNoGTO => 1);
        } else {
            my $geo = GEO->CreateFromGto($gtoFile, detail => 0);
            if (! $geo->good_seed) {
                $stats->Add(genomeBadSeed => 1);
            } else {
                # Get the FASTA file.
                my $fastaFile = $fastaFiles{$gName};
                if ($fastaFile) {
                    $stats->Add(genomeFasta => 1);
                } else {
                    $fastaFile = '';
                }
                P3Utils::print_cols([$genomeID, $genomeName, $fastaFile, $gtoFile, $pStatus]);
            }
        }
    }
    $count++;
    if ($debug && $count % 1000 == 0) {
        print STDERR "$count genomes processed.\n";
    }
}
# All done.  Show the statistics.
print STDERR "All done.\n" . $stats->Show() if $debug;

## Get all the genome IDs currently in PATRIC.  We map them to 'public' or 'private'.
sub get_all_genomes {
    my ($p3) = @_;
    print STDERR "Retrieving IDs of genomes already in PATRIC.\n" if $debug;
    my $results = P3Utils::get_data($p3, genome => [['ne', 'genome_id', 0]], ['genome_id', 'public']);
    my %retVal = map { $_->[0] => ($_->[1] ? 'public' : 'private') } @$results;
    print STDERR scalar(keys %retVal) . " genome IDs found.\n" if $debug;
    return \%retVal;
}