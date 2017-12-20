=head1 Preliminary Check for Good Genomes

    p3-good-check.pl [options] inDir outDir

This script performs a prima facie good-genome check on the PATRIC genomes. It accepts as input a list of previously-determined
good and bad genome IDs. It first checks for a bad seed protein, then runs through a SciKit check on everything remaining.
It will output a list of newly-discovered bad genomes and a list of genomes that need to be run through CheckM.

=head2 Parameters

The positional parameters are the names of the input and output directories.

The input directory must contain a C<good.patric.tbl> file containing the known good-genome IDs, and a C<bad.patric.tbl> file containing the known
bad-genome IDs. New versions of these files will be put into the output directory and the list of genomes to check will be written to C<check.tbl>.

The command-line options are as follows:

=over 4

=item genomes

If specified, the name of a tab-delimited file containing (0) the ID and (1) the name of each genome to check. The default is to check all public genomes.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use IO::Handle;
use File::Copy::Recursive;
use SeedUtils;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir',
        ['genomes=s', 'check only genomes in the specified file']);
# Get the working directory.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    die "Invalid or missing output directory $outDir.";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# To start, we get the known-genome information.
my $goodH = read_ids("$inDir/good.patric.tbl");
my $badH = read_ids("$inDir/bad.patric.tbl");
# Now get all of the genome IDs and names.
my $genomeList = [];
if ($opt->genomes) {
    # Here we have a file of IDs and names.
    print "Reading genomes from file.\n";
    open(my $ih, '<', $opt->genomes) || die "Could not open genome file: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\S+)\t(.+)/) {
            push @$genomeList, [$1, $2];
        }
    }
} else {
    print "Reading genomes from PATRIC.\n";
    $genomeList = P3Utils::get_data($p3, genome => [['eq', 'public', 1]], ['genome_id', 'genome_name']);
}
my $total = scalar(@$genomeList);
print "$total genomes to check.\n";
# Sort the output.
$genomeList = [ sort { $a->[0] cmp $b->[0] } @$genomeList ];
# Create the temp directory for the SciKit tool.
my $pDir = "$outDir/Temp";
if (! -d $pDir) {
    print "Creating temp directories for SciKit.\n";
    File::Copy::Recursive::pathmk($pDir);
    File::Copy::Recursive::pathmk("$pDir/SciKit");
}
# Open the output files.
my $gh = IO::File->new(">$outDir/good.patric.tbl") || die "Could not open output good.patric.tbl: $!";
my $bh = IO::File->new(">$outDir/bad.patric.tbl") || die "Could not open output bad.patric.tbl: $!";
my $ch = IO::File->new(">$outDir/check.tbl") || die "Could not open output check.tbl: $!";
# This list will accumulate a batch of genomes of unknown provenance.
my @queue;
# Run through the PATRIC genomes.
for my $genomeEntry (@$genomeList) {
    my ($id, $name) = @$genomeEntry;
    # Check the hashes first.
    if ($goodH->{$id}) {
        record($genomeEntry, $gh);
    } elsif ($badH->{$id}) {
        record($genomeEntry, $bh);
    } else {
        # Here we need to queue up the genome for further processing.
        push @queue, $genomeEntry;
    }
}
# Process the residual.
print "Processing residual.  " . scalar(@queue) . " left to check.\n";
my $n = scalar @queue - 1;
for (my $i = 0; $i <= $n; $i += 100) {
    my $j = $i + 99; $j = $n if $j > $n;
    my $subset = [@queue[$i .. $j]];
    process_batch($subset);
}
print "All done.\n";


# Record a genome entry in a file.
sub record {
    my ($gEntry, $fh) = @_;
    my $line = join("\t", @$gEntry) . "\n";
    $fh->printflush($line);
}


# Process a batch of unknown-quality genomes.
sub process_batch {
    my ($queue) = @_;
    # First, check the seed proteins.
    my @genomes = map { $_->[0] } @$queue;
    print "Checking seed proteins for " . scalar(@genomes) . " genomes.\n";
    my $features = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', 'Phenylalanyl tRNA-synthetase alpha chain']],
        ['genome_id', 'aa_sequence'], \@genomes, 'genome_id');
    print "Analyzing seed protein lengths.\n";
    my %lengths = map { $_ => [] } @genomes;
    for my $feature (@$features) {
        my ($genome, $seq) = @$feature;
        push @{$lengths{$genome}}, length($seq);
    }
    print "Processing individual genomes.\n";
    for my $genomeData (@$queue) {
        my ($genome, $name) = @$genomeData;
        my $lengthList = $lengths{$genome};
        if (scalar @$lengthList != 1) {
            record($genomeData, $bh);
            print "$genome $name has a bad protein count.\n";
        } else {
            my $aaLen = $lengthList->[0];
            if ($aaLen < 209 || $aaLen > 405) {
                record($genomeData, $bh);
                print "$genome $name has a bad protein length.\n";
            } else {
                # We need to do a quality check here. Get the GTO and write its
                # JSON to disk.
                my $start = time;
                print "Retrieving GTO for $genome $name.\n";
                my $gto = $p3->gto_of($genome);
                $gto->destroy_to_file("$pDir/bin.gto");
                undef $gto;
                # Clean up past working files from the checkers.
                print "Cleaning work directories.\n";
                File::Copy::Recursive::pathempty("$pDir/SciKit") || die "Could not clean SciKit working directory: $!";
                my $cmd = "gto_consistency $pDir/bin.gto $pDir/SciKit $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $FIG_Config::global/roles.to.use";
                SeedUtils::run($cmd);
                my $score = 0;
                if (! open(my $ih, '<', "$pDir/SciKit/evaluate.log")) {
                    print "WARNING: Cannot open output from Scikit: $!\n";
                } else {
                    while (! eof $ih) {
                        my $line = <$ih>;
                        if ($line =~ /^Fine_Consistency=\s+(.+)%/) {
                            $score = $1;
                        }
                    }
                }
                print "SciKit fine score is $score. " . (time - $start) . " seconds to check.\n";
                if ($score < 85) {
                    print "$genome $name rejected by SciKit.\n";
                    record($genomeData, $bh);
                } else {
                    print "$genome $name queued for CheckM.\n";
                    record($genomeData, $ch);
                }
            }
        }
    }
}


# Read genome IDs from a file.
sub read_ids {
    my ($fileName) = @_;
    my %retVal;
    my $count = 0;
    open(my $ih, '<', $fileName) || die "Could not open $fileName: $!";
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /(\d+\.\d+)/) {
            $retVal{$1} = 1;
            $count++;
        }
    }
    print "$count genomes read from $fileName.\n";
    return \%retVal;
}
