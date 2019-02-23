=head1 Check for Good PATRIC Genomes

    p3-good-check.pl [options] inDir outDir

This script performs a good-genome check on the PATRIC genomes. It accepts as input a list of previously-determined
good and bad genome IDs. It first checks for a bad seed protein, then performs an EvalG check, then runs through a
SciKit check on everything remaining.

It will output a list of newly-discovered good and bad genomes.

=head2 Parameters

The positional parameters are the names of the input and output directories.

The input directory must contain a C<good.patric.tbl> file containing the known good-genome IDs and a C<bad.patric.tbl> file
containing the known bad-genome IDs. New versions of these files will be put into the output directory. If the input directory
is C<null>, then the output directory will be created fresh.

The command-line options are as follows:

=over 4

=item genomes

If specified, the name of a tab-delimited file containing the ID of each genome to check in the first column.
The default is to check all public genomes.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use IO::Handle;
use File::Copy::Recursive;
use SeedUtils;
use EvalCom::Tax;
use RoleParse;
use Stats;
use GPUtils;

$| = 1;
my $stats = Stats->new();
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir',
        ['genomes=s', 'check only genomes in the specified file']);
# Get the working directory.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif ($inDir ne 'null' && ! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory: $!";
}
# Get access to PATRIC.
print "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# To start, we get the known-genome information.
my ($goodH, $badH) = ({}, {});
if ($inDir ne 'null') {
    $goodH = read_ids("$inDir/good.patric.tbl");
    $badH = read_ids("$inDir/bad.patric.tbl");
}
# Now get all of the genome IDs and names.
my $genomeList = [];
if ($opt->genomes) {
    # Here we have a file of IDs and names.
    print "Reading genomes from file.\n";
    open(my $ih, '<', $opt->genomes) || die "Could not open genome file: $!";
    my @gids;
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^(\d+\.\d+)/) {
            push @gids, $1;
        }
    }
    print "Getting " . scalar(@gids) . " genomes from PATRIC.\n";
    $genomeList = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name', 'kingdom'], \@gids, 'genome_id');
} else {
    print "Reading genomes from PATRIC.\n";
    $genomeList = P3Utils::get_data($p3, genome => [['eq', 'public', 1]], ['genome_id', 'genome_name', 'kingdom']);
}
my $total = scalar(@$genomeList);
print "$total genomes to check.\n";
# Sort the output.
$genomeList = [ sort { $a->[0] cmp $b->[0] } @$genomeList ];
# Get the completeness checker.
print "Initializing GTO checker.\n";
my $checkG = EvalCom::Tax->new("$FIG_Config::global/CheckG", logH => \*STDOUT, stats => $stats);
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
# This list will accumulate a batch of genomes of unknown provenance.
my @queue;
# Run through the PATRIC genomes.
for my $genomeEntry (@$genomeList) {
    my ($id, $name, $type) = @$genomeEntry;
    # Only keep proks.
    $stats->Add($type => 1);
    if ($type eq 'Archaea' || $type eq 'Bacteria') {
        # Check the hashes first.
        if ($goodH->{$id}) {
            record($genomeEntry, $gh);
            $stats->Add(oldGood => 1);
        } elsif ($badH->{$id}) {
            record($genomeEntry, $bh);
            $stats->Add(oldBad => 1);
        } else {
            # Here we need to queue up the genome for further processing.
            push @queue, $genomeEntry;
            $stats->Add(newGenome => 1);
        }
    }
}
# Process the residual.
print "Processing queue.  " . scalar(@queue) . " to check.\n";
my $n = scalar @queue - 1;
for (my $i = 0; $i <= $n; $i += 100) {
    my $j = $i + 99; $j = $n if $j > $n;
    my $subset = [@queue[$i .. $j]];
    process_batch($subset);
    $stats->Add(batchesCompleted => 1);
    my $remaining = $n - $j;
    print "$remaining genomes left to check.\n";
}
print "All done.\n" . $stats->Show();


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
        ['genome_id', 'product', 'aa_sequence'], \@genomes, 'genome_id');
    print "Analyzing seed protein lengths.\n";
    my %lengths = map { $_ => [] } @genomes;
    for my $feature (@$features) {
        my ($genome, $function, $seq) = @$feature;
        my $checksum = RoleParse::Checksum($function);
        if ($checksum eq 'WCzieTC/aZ6262l19bwqgw') {
            push @{$lengths{$genome}}, length($seq);
        }
    }
    print "Processing individual genomes.\n";
    for my $genomeData (@$queue) {
        my ($genome, $name) = @$genomeData;
        $stats->Add(genomeIn => 1);
        my $lengthList = $lengths{$genome};
        if (scalar @$lengthList != 1) {
            record($genomeData, $bh);
            print "$genome $name has a bad seed protein count.\n";
            $stats->Add(badSeedCount => 1);
        } else {
            my $start = time;
            my $bad = 0;
            my $aaLen = $lengthList->[0];
            if ($aaLen < 209 || $aaLen > 652) {
                $bad = 1;
                print "$genome $name has a bad seed protein length.\n";
                $stats->Add(badSeedLength => 1);
            } else {
                # We need to do a quality check here. Get the GTO. We have to do a full good-seed
                # check (the above was just an easy filter, since we don't know domain yet), a
                # completeness check, and if it passes all that, write its JSON to disk for a
                # consistency check.
                print "Retrieving GTO for $genome $name.\n";
                my $gto = $p3->gto_of($genome);
                print "Checking completeness.\n";
                my $evalH = $checkG->Check($gto);
                my ($complete, $contam) = ($evalH->{complete}, $evalH->{contam});
                if (! defined $complete) {
                    print "$genome $name is not prokaryotic.\n";
                    $stats->Add(notProk => 1);
                    $bad = 1;
                } else {
                    print "$genome $name has completeness $complete and $contam contamination.\n";
                    if ($complete < 80) {
                        $bad = 1;
                        $stats->Add(incomplete => 1);
                        print "$genome $name is incomplete.\n";
                    } elsif ($contam > 10) {
                        $bad = 1;
                        $stats->Add(contaminated => 1);
                        print "$genome $name is contaminated.\n";
                    } elsif (! GPUtils::good_seed($gto)) {
                        $bad = 1;
                        $stats->Add(badSeedLength => 1);
                        print "$genome $name has a bad seed protein length.\n";
                    } else {
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
                        print "SciKit fine score is $score.\n";
                        if ($score < 87) {
                            print "$genome $name rejected by SciKit.\n";
                            $stats->Add(inconsistent => 1);
                            $bad = 1;
                        }
                    }
                }
            }
            if ($bad) {
                record($genomeData, $bh);
            } else {
                record($genomeData, $gh);
            }
            print "$genome checked in " . (time - $start) . " seconds.\n";
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
