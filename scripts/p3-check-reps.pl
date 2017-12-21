
=head1 Check For Representative Genomes

    p3-check-reps.pl [options] inDir outDir

This script looks at an input list of PATRIC genome IDs and determines the representative genome for each one.

The script operates on a representative server directory. The directory contains four files of paramount interest.

=over 4

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item rep_db.tbl

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID, (1) the ID of the genome's
representative, and (2) the similarity number.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

=head2 Parameters

The positional parameters are the input directory and the output directory. The input directory must contain the above three files.

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain PATRIC genome IDs in the key column.

Additional command-line options are those given in L<P3Utils/col_options> (to select the input column) plus the following.

=over 4

=item checkOnly

If specified, no new representative genomes will be computed. The genomes that were unrepresented will be output to the file
C<outliers.tbl> in the output directory.

=item filter

If specified, the outlier genomes will be filtered so that only good ones are checked for representatives.

=item middle

If specified, a similarity score indicating a middle similarity distance. The file C<middle.tbl> will be written containing
the genomes that are less similar than the main score but more similar than this score.

=back

=head2 Output Files

The standard output will contain a progress report.

The output directory will contain an updated version of the input directory with as many genomes represented as possible.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir', P3Utils::col_options(), P3Utils::ih_options(),
        ['checkOnly', 'check for representation only-- do not add new representative genomes'],
        ['filter', 'perform a good-genome analysis on the unrepresented genomes'],
        ['middle=i', 'alternate similarity score for use in determining mid-distance genomes']
        );
# Create the statistics object.
my $stats = Stats->new();
# Get access to PATRIC.
print "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# Verify the parameters.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! -f "$inDir/6.1.1.20.fasta" || ! -f "$inDir/complete.genomes") {
    die "$inDir does not appear to contain a representative-genomes database.";
} else {
    print "Input database in $inDir.\n";
}
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} else {
    print "Output directory is $outDir.\n";
}
# Create the PATRIC filter and column clauses for genome queries.
my @filter = (['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']);
my @cols = qw(genome_id genome_name aa_sequence);
# Create the database from the input directory.
print "Creating database from $inDir.\n";
my $repDB = RepGenomeDb->new_from_dir($inDir, verbose => 1);
# Save the parameters.
my $K = $repDB->K();
my $minScore = $repDB->score();
print "Kmer size is $K and minimum similarity is $minScore.\n";
# Compute the middle distance.
my $midScore = $opt->middle || $minScore;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# This will be a list of RepGenome objects for the un-represented genomes.
my @outliers;
# This will be a list of RepGenome objects for un-represented genomes at the middle distance.
my @middle;
# This will be a hash of bad genome IDs.
my %bad;
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    my @genomes = map { $_->[0] } @$couplets;
    # We need to isolate the genomes for which we don't already have a representative.
    my @lost;
    for my $genome (@genomes) {
        $stats->Add(genomeIn => 1);
        my ($repID, $score) = $repDB->check_rep($genome);
        if (! $repID) {
            push @lost, $genome;
            $stats->Add(genomeNotFound => 1);
        } else {
            $stats->Add(genomeFound => 1);
        }
    }
    if (@lost) {
        print scalar(@lost) . " genomes in this batch are not yet in the database.\n";
        # Ask PATRIC for the names and identifying proteins of the un-represented genomes.
        my $resultList = P3Utils::get_data_keyed($p3, 'feature', \@filter, \@cols, \@lost, 'genome_id');
        # The resultList entries are in the form [$genome, $name, $prot]. Get the longest
        # protein for each genome. We also track obviously bad genomes.
        my (%results);
        for my $result (@$resultList) {
            my ($genome, $name, $prot) = @$result;
            # Check the protein.
            if (! $prot) {
                print "WARNING: $genome $name has no identifying protein.\n";
                $stats->Add(genomeNoProt => 1);
                $stats->Add(badGenome => 1);
                $bad{$genome} = 1;
            } else {
                # Add the protein length to the result array.
                my $protLen = length $prot;
                push @$result, $protLen;
                if (! exists $results{$genome}) {
                    $results{$genome} = $result;
                } else {
                    $stats->Add(redundantProt => 1);
                    if (! $bad{$genome}) {
                        print "WARNING: $genome $name has a redundant identifying protein.\n";
                        $bad{$genome} = 1;
                        $stats->Add(badGenome => 1);
                    }
                    if ($protLen > $results{$genome}[3]) {
                        # It's a better protein, so keep it.
                        $results{$genome} = $result;
                    }
                }
            }
        }
        # Now loop through the genomes, checking for representatives.
        for my $genome (keys %results) {
            print "Checking $genome.\n";
            my (undef, $name, $prot) = @{$results{$genome}};
            my ($repID, $score) = $repDB->find_rep($prot);
            if ($score >= $minScore) {
                print "$genome $name assigned to $repID with similarity $score.\n";
                $repDB->Connect($repID, $genome, $score);
                $stats->Add(genomeConnected => 1);
            } else {
                my $repGenome = RepGenome->new($genome, name => $name, prot => $prot, K => $K);
                if ($score >= $midScore) {
                    # Here we are at the middle distance.
                    print "$genome $name is at the middle distance from $repID with similarity $score.\n";
                    push @middle, $repGenome;
                    $stats->Add(genomeMiddle => 1);
                } else {
                    print "$genome $name is an outlier. Best score was $score.\n";
                    $stats->Add(genomeOutlier => 1);
                }
                push @outliers, $repGenome;
            }
        }
    }
}
# Checkpoint our results.
print "Saving results of initial run.\n";
$repDB->Save($outDir);
if ($opt->middle) {
    dump_outliers("$outDir/middle.tbl", \@middle);
}
# Now we need to process the outliers.
my $nOutliers = scalar @outliers;
my $mCount = scalar @middle;
my $nCount = $nOutliers - $mCount;
print "$nOutliers un-represented genomes found.\n";
if ($opt->middle) {
    print "$mCount at the middle distance, $nCount true outliers.\n";
}
if ($opt->filter) {
    print "Filtering for bad outliers.\n";
    my $pDir = "$outDir/Temp";
    if (! -d $pDir) {
        File::Copy::Recursive::pathmk($pDir);
    }
    my @good;
    for my $outlier (@outliers) {
        my $id = $outlier->id();
        my $name = $outlier->name();
        if ($bad{$id}) {
            print "$id $name has too many seed proteins.\n";
        } else {
            # Check the protein length.
            my $aaLen = length $outlier->prot();
            if ($aaLen < 209 || $aaLen > 405) {
                print "$id $name has a bad protein length.\n";
            } else {
                # We need to do a quality check here. Get the GTO and write its FASTA
                # and JSON to disk.
                print "Retrieving GTO for $id $name.\n";
                my $gto = $p3->gto_of($id);
                $gto->write_contigs_to_file("$pDir/bin.fa");
                $gto->destroy_to_file("$pDir/bin.gto");
                undef $gto;
                # Clean up past working files from the checkers.
                print "Cleaning work directories.\n";
                if (-d "$pDir/Eval/SciKit") {
                    File::Copy::Recursive::pathempty("$pDir/Eval/SciKit") || die "Could not clean SciKit working directory: $!";
                } else {
                    File::Copy::Recursive::pathmk("$pDir/Eval/SciKit") || die "Could not create SciKit working directory: $!";
                }
                if (-d "$pDir/Eval/CheckM") {
                    File::Copy::Recursive::pathempty("$pDir/Eval/CheckM") || die "Could not clean CheckM working directory: $!";
                } else {
                    File::Copy::Recursive::pathmk("$pDir/Eval/CheckM") || die "Could not create CheckM working directory: $!";
                }
                # We do the SciKit check first, because it's faster.
                my $cmd = "gto_consistency $pDir/bin.gto $pDir/Eval/SciKit $FIG_Config::global/FunctionPredictors $FIG_Config::global/roles.in.subsystems $FIG_Config::global/roles.to.use";
                SeedUtils::run($cmd);
                $stats->Add(sciKitRun => 1);
                my $score = 0;
                if (! open(my $fh, '<', "$pDir/Eval/SciKit/evaluate.log")) {
                    print "WARNING: Cannot open output from Scikit: $!\n";
                } else {
                    while (! eof $fh) {
                        my $line = <$fh>;
                        if ($line =~ /^Fine_Consistency=\s+(.+)%/) {
                            $score = $1;
                        }
                    }
                }
                print "SciKit fine score is $score.\n";
                if ($score >= 85) {
                    # We have passed SciKit. Now try CheckM.
                    $cmd = "checkm lineage_wf --tmpdir $FIG_Config::temp -x fa --file $pDir/checkm.log $pDir $pDir/Eval/CheckM";
                    SeedUtils::run($cmd);
                    $stats->Add(checkMrun => 1);
                    my ($checkMscore, $checkMcontam) = (0, 999);
                    if (! open(my $fh, '<', "$pDir/checkm.log")) {
                        print "WARNING: Cannot open output from CheckM: $!";
                    } else {
                        while (! eof $fh) {
                            my $line = <$fh>;
                            if ($line =~ /^\s+bin\s+/) {
                                my @cols = split /\s+/, $line;
                                $checkMscore = $cols[13];
                                $checkMcontam = $cols[14];
                            }
                        }
                        close $fh;
                    }
                    print "CheckM score is $checkMscore complete, $checkMcontam contamination.\n";
                    if ($checkMscore >= 80 && $checkMcontam <= 10) {
                        push @good, $outlier;
                        $stats->Add(goodOutlier => 1);
                    }
                }
            }
        }
    }
    # End of outlier loop.
    @outliers = @good;
    $nOutliers = scalar @outliers;
    print "$nOutliers outlier genomes are good.\n";
}
if ($opt->checkonly) {
    # We don't want to update the representatives, we just want a list of outliers.
    open(my $oh, '>', "$outDir/outliers.tbl") || die "Could not open outliers file: $!";
    print "Saving list of outliers.\n";
    dump_outliers("$outDir/outliers.tbl", \@outliers);
} else {
    # Loop through the outliers, computing a cross-reference. The cross-reference is a two-dimensional table
    # of similarity scores, indexed by the genome's position in the outlier list.
    my $xref = [];
    my ($i, $j);
    for ($i = 0; $i < $nOutliers; $i++) {
        my $repGenome = $outliers[$i];
        print "Cross-referencing " . $repGenome->id() . "\n";
        # Copy the scores already computed.
        for ($j = 0; $j < $i; $j++) {
            $xref->[$i][$j] = $xref->[$j][$i];
        }
        # Compute the remaining scores.
        for ($j = $i; $j < $nOutliers; $j++) {
            my $prot = $outliers[$j]->prot();
            $xref->[$i][$j] = $repGenome->check_genome($prot);
        }
    }
    # Now we process the scores until we run out of outliers.
    # We checkpoint each time this counter drops to zero.
    my $lastCheck = 1000;
    # We stop everything when this switch is tripped.
    my $done = 0;
    # Here is the loop.
    while (@outliers && ! $done) {
        print "Searching for most popular genome among $nOutliers outliers.\n";
        # Find the genome with the most similarities.
        my ($b, $best) = (0, count($xref->[0], $minScore));
        for ($i = 1; $i < $nOutliers; $i++) {
            my $count = count($xref->[$i], $minScore);
            if ($count > $best) {
                $b = $i; $best = $count;
            }
        }
        my $bestGenome = $outliers[$b];
        my $bestID = $bestGenome->id();
        my $bestName = $bestGenome->name();
        if ($best <= 1) {
            print "End of useful genomes.\n";
            $done = 1;
        } else {
            print "Best genome is $bestID $bestName with $best neighbors.\n";
            # This list will contain the indices of the outliers we are keeping; that is, the
            # ones not similar to the best one.
            my @keepers;
            my $connected = 0;
            # Loop through the outliers, connecting the similar genomes to the best one.
            for ($i = 0; $i < $nOutliers; $i++) {
                my $score = $xref->[$b][$i];
                if ($score >= $minScore) {
                    $bestGenome->AddGenome($outliers[$i]->id(), $score);
                    $stats->Add(genomeConnected => 1);
                    $connected++;
                } else {
                    push @keepers, $i;
                }
            }
            print "$connected genomes represented by $bestID.\n";
            # Add the best genome to the database.
            $repDB->AddRepObject($bestGenome);
            $stats->Add(newRepAdded => 1);
            # Form the new outliers list.
            @outliers = map { $outliers[$_] } @keepers;
            $nOutliers = scalar @outliers;
            # Check for a checkpoint.
            $lastCheck -= $connected;
            if ($lastCheck <= 0) {
                print "Saving results to $outDir.\n";
                $repDB->Save($outDir);
                $lastCheck = 1000;
            }
            # Form the new xref. We map the old indices to the new ones, which eliminates the
            # genomes we connected above.
            my @newXref;
            print "Reorganizing the remaining $nOutliers genomes.\n";
            for ($i = 0; $i < $nOutliers; $i++) {
                for ($j = 0; $j < $nOutliers; $j++) {
                    $newXref[$i][$j] = $xref->[$keepers[$i]][$keepers[$j]];
                }
            }
            $xref = \@newXref;
        }
    }
    # All the outliers have been added to the database.
    print "Saving to $outDir.\n";
    $repDB->Save($outDir);
    if (@outliers) {
        print "Printing useless outliers.\n";
        dump_outliers("$outDir/useless.tbl", \@outliers);
        $stats->Add(unprocessedGenome => scalar @outliers);
    }
}
print "All done.\n" . $stats->Show();


# Write out the list of remaining outliers.
sub dump_outliers {
    my ($fileName, $outliers) = @_;
    open (my $oh, '>', $fileName) || die "Could not open outlier file: $!";
    print $oh "id\tname\tprot\n";
    for my $outlier (@$outliers) {
        print $oh join("\t", $outlier->id(), $outlier->name(), $outlier->prot()) . "\n";
    }
}

# Count the good similarities in a list of scores.
sub count {
    my ($list, $minScore) = @_;
    my $retVal = 0;
    for my $score (@$list) {
        if ($score >= $minScore) {
            $retVal++;
        }
    }
    return $retVal;
}



