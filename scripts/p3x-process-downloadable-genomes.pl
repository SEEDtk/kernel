=head1 Submit Downloadable Genomes to PATRIC

    p3x-process-downloadable-genomes.pl [options] workDir inFile folder

This script will read a file of genome data, including FASTA file URLs and taxonomic information.  The FASTA files will be
downloaded into the specified directory and then submitted to PATRIC RAST.  After the jobs complete, the results will
be interrogated for quality information.

Status messages will be written to the standard error output.  This should be redirected separately from the standard output.

=head2 Parameters

The positional parameters are the name of the working directory, the name of the input file, and the name of the output folder
in the user's PATRIC workspace.  The user must be logged into PATRIC using L<p3-login.pl>.

The output folder (if specified) should be a fully-qualified name (e.g. C</user@patricbrc.org/homes/DownloadedGenomes>).  The
genome will be put in a sub-folder with the same name as the genome's name, computed from the species name and the label.  If
a genome with that name already exists in the output folder, an internal server error will occur and the entire script will fail.
Thus, this capability is best used sparingly.

Additional command-line options are L<P3Utils/col_options>.  The key column in this case is the one containing the FASTA file
URLs.  The batch size is the number of genomes to submit at one time.  The script will wait for all of them to finish before
submitting more, as a throttling measure.

The command-line options in L<Shrub/script_options> may be used to specify the L<Shrub> database, which is used to convert
species names to taxonomic grouping IDs.

The following additional command-line options are supported.

=over 4

=item label

The position (1-based) or name of the input column containing the label to give to the new genome.  The default is C<1>,
indicating the first column.

=item retain

A comma-delimited list of the positions (1-based) or names of input columns to retain in the output.  The default is only the
label column retained.  The label is not retained automatically if this option is non-empty.

=item species

The position (1-based) or name of the column containing the species name.  The default is none, indicating the species must be
estimated from the contigs.

=item skip

The number of genomes to skip before starting.  This allows you to resume after an error, or to split a large file into multiple
runs. The default is C<0>

=item process

The number of genomes to process.  The default value of C<0> indicates processing should continue to the end of the file.

=item gtos

If specified, the L<GenomeTypeObject> files for the genomes created will be written to the working directory.  Otherwise,
they will be discarded.

=item clear

If specified, and the working directory exists, it will be erased before we start.

=item force

If specified, then when downloading, the file will always be downloaded.  The default behavior is to skip downloading if the
output file already exists.

=item noindex

Suppress indexing of the genome in the PATRIC database.  Use this to keep the genomes from showing up on the PATRIC web site or via the
CLI.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Shrub;
use TaxCheck;
use RASTlib;
use File::Copy::Recursive;
use LWP::Simple;
use GEO;

# Get the command-line options.
my $opt = P3Utils::script_opts('workDir folder', P3Utils::col_options(), Shrub::script_options(),
        ['label=s', 'input column containing the genome label', { default => 1 }],
        ['retain=s', 'comma-delimited list of columns to retain'],
        ['species=s', 'input column containing the genome species name'],
        ['skip=i', 'number of input genomes to skip', { default => 0 }],
        ['process=i', 'number of input genomes to process (0 for all)', { default => 0 }],
        ['gtos', 'save genome GTO files to disk'],
        ['clear', 'erase the working directory before starting'],
        ['force', 'always download FASTA files'],
        ['noindex', 'do not put the genomes into the PATRIC index']
        );
# Get the main parameters.
my ($workDir, $inFile, $folder) = @ARGV;
if (! $workDir) {
    die "No working directory specified.";
} elsif (! -d $workDir) {
    print STDERR "Creating work directory $workDir.\n";
    File::Copy::Recursive::pathmk($workDir) || die "Could not create $workDir: $!";
} elsif ($opt->clear) {
    print STDERR "Erasing work directory $workDir.\n";
    File::Copy::Recursive::pathempty($workDir) || die "Could not erase $workDir: $!";
}
my $ih;
if (! $inFile) {
    die "No input file specified.";
} elsif (! -s $inFile) {
    die "Input file $inFile not found or empty.";
} else {
    open($ih, '<', $inFile) || die "Could not open $inFile: $!";
}
if (! $folder) {
    print STDERR "Default folder QuickData selected for workspace output.\n";
} else {
    print STDERR "Workspace output will be put in $folder.\n";
}
# Get access to PATRIC.
print STDERR "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# Connect to the Shrub database.
print STDERR "Connecting to Shrub.\n";
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
print STDERR "Parsing input headers.\n";
# Read the incoming headers.
my ($inHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Locate the columns to retain.
my $label = $opt->label;
my $retainCols = ($opt->retain // [$label]);
(undef, $retainCols) = P3Utils::find_headers($inHeaders, input => split(/,/, $retainCols));
# Find the label column.
my $labelCol = P3Utils::find_column($label, $inHeaders);
print STDERR "Label is in column " . ($labelCol + 1) . "\n";
# If there is a species column, locate that, too.
my $speciesCol;
if ($opt->species) {
    $speciesCol = P3Utils::find_column($opt->species, $inHeaders);
    print STDERR "Species information will be taken from input file.\n";
} else {
    print STDERR "Species information will be estimated from FASTA files.\n";
}
# Create the taxonomy checker.  Note we default on everything.
my $taxCheck = TaxCheck->new($p3, debug => 1);
# Create the authorization header.  Note we are assuming the user is logged in, because we are not
# passing in credentials.
my $header = RASTlib::auth_header();
# Do any skipping we need to do.  Note the header has already been read.
my $skip = $opt->skip;
if ($skip) {
    print STDERR "Skipping $skip lines in file.\n";
}
while ($skip > 0) {
    my $line = <$ih>;
    $skip--;
}
# This will accumulate the jobs we could not process.
my @failures;
# Finally, initialize the output headers.
if (! $opt->nohead) {
    my @outHeaders;
    for my $retainCol (@$retainCols) {
        push @outHeaders, $inHeaders->[$retainCol];
    }
    push @outHeaders, 'PATRIC ID', 'Scientific Name', 'EvalG Completeness', 'EvalG Contamination', 'Coarse Consistency', 'Fine Consistency', 'Good?';
    P3Utils::print_cols(\@outHeaders);
}
# Compute the process number.
my $process = $opt->process;
if ($process == 0) {
    print STDERR "Processing all records.\n";
} else {
    print STDERR "Stopping at $process records.\n";
}
# Loop through the input.
my $lNum = 0;
my $done;
my $processed = 0;
my $goodFound = 0;
my $start = time;
while (! $done) {
    # Get the first batch of records.
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    print STDERR scalar(@$couplets) . " lines in input batch.\n";
    # We will process each of the genomes individually.  At the end, we run through them asking for status.
    # The following hashes are keyed by genome label.  This maps the label to the proposed output line of
    # retained input columns.  Later, the output data [genomeID, name, completeness, contamination, coarseConsistency, fineConsistency]
    # will be appended to each line.
    my %gRetain;
    # This maps the label to the RAST job ID.
    my %gJobs;
    # Loop through the genomes, submitting.
    for my $couplet (@$couplets) {
        if (! $done) {
            my ($url, $line) = @$couplet;
            $lNum++;
            # Extract the label.
            my $label = $line->[$labelCol];
            print STDERR "Processing line $lNum: $label.\n";
            # Get the retained columns.
            my @retainers = P3Utils::get_cols($line, $retainCols);
            # The next step is to download the FASTA file.
            my $fastaFile = "$workDir/$label.fa";
            my $fileOK;
            if (-s $fastaFile && ! $opt->force) {
                print STDERR "$fastaFile already downloaded.\n";
                $fileOK = 1;
            } else {
                print STDERR "Downloading $url into $fastaFile.\n";
                my $rc = getstore($url, $fastaFile);
                if (is_success($rc)) {
                    $fileOK = 1;
                } else {
                    print STDERR "DOWNLOAD ERROR for $label: status = $rc.\n";
                    push @failures, $label;
                }
            }
            if ($fileOK) {
                # Compute the taxonomic information.
                my ($domain, $species, $speciesID);
                if (defined $speciesCol) {
                    # Here we get the species name from the input.
                    $species = $line->[$speciesCol];
                }
                # If the species column is nonblank, use it.
                if ($species && $species !~ /^n\/a/i) {
                    # First we get the ID.
                    ($speciesID) = $shrub->GetFlat('TaxonomicGrouping', 'TaxonomicGrouping(scientific-name) = ?', [$species], 'id');
                    if (defined $speciesID) {
                        # Now the domain.
                        $domain = $shrub->domain_of($speciesID);
                    }
                }
                # Here we have to estimate the species.
                if (! defined $speciesID) {
                    my ($newID, $newName, $lineage) = $taxCheck->Compute($fastaFile, 'species');
                    if (! $newID) {
                        print STDERR "COULD NOT COMPUTE SPECIES for $fastaFile.\n";
                        push @failures, $label;
                    } else {
                        $speciesID = $newID;
                        $species = $newName;
                        $domain = $lineage->[0];
                    }
                }
                if (defined $speciesID) {
                    # Convert the domain to a domain code.
                    $domain = uc substr($domain, 0, 1);
                    # We are about to submit.  Save the retained columns.
                    $gRetain{$label} = \@retainers;
                    # Submit the job to RAST.
                    print STDERR "Submitting $fastaFile using $species and domain $domain.\n";
                    my $contigs = RASTlib::read_fasta($fastaFile);
                    $gJobs{$label} = RASTlib::Submit($contigs, $speciesID, "$species $label", domain => $domain,
                            path => $folder, header => $header, noIndex => $opt->noindex);
                    print STDERR "Job ID is $gJobs{$label}.\n";
                }
            }
            if ($process == 1) {
                $done = 1;
            } else {
                $process--;
            }
        }
    }
    # Now all the jobs have been submitted.  The next step is to wait for them to finish.
    my @incomplete = keys %gJobs;
    my $count = scalar @incomplete;
    while ($count) {
        print STDERR "Waiting...\n";
        sleep 10;
        my @remain;
        print STDERR "Checking status of $count jobs.\n";
        for my $job (@incomplete) {
            my $jobID = $gJobs{$job};
            if (! RASTlib::check($jobID, $header)) {
                push @remain, $job;
            } else {
                print STDERR "$job has completed.\n";
                my $gto = RASTlib::retrieve($jobID, $header);
                # Get the genome ID and name.
                my $genomeID = $gto->{id};
                my $genomeName = $gto->{scientific_name};
                # Check the GTOS option.
                if ($opt->gtos) {
                    # Save the GTOs to disk.
                    open(my $oh, '>', "$workDir/$genomeID.gto") || die "Could not open GTO file for $genomeID: $!";
                    print $oh $gto;
                    close $oh;
                    print STDERR "Genome saved to $genomeID.gto.\n";
                }
                my $quality = $gto->{quality};
                my $fine = $quality->{fine_consistency};
                my $coarse = $quality->{coarse_consistency};
                my $complete = $quality->{completeness};
                my $contam = $quality->{contamination};
                my $good = ((GEO::completeX($complete) && GEO::contamX($contam) && GEO::consistX($fine)) ? 'Y' : '');
                $goodFound++ if $good;
                push @{$gRetain{$job}}, $genomeID, $genomeName, $complete, $contam, $coarse, $fine, $good;
                $count--;
                $processed++;
            }
        }
        @incomplete = @remain;
    }
    # Write the output rows.
    print STDERR "Writing to output.\n";
    for my $job (sort keys %gJobs) {
        P3Utils::print_cols($gRetain{$job});
    }
    # Stop if there is no more input.
    if (eof $ih) {
        $done = 1;
    }
}
my $failCount = scalar @failures;
if ($failCount) {
    print STDERR "Failed to run $failCount jobs.\n";
    for my $failure (@failures) {
        print STDERR "     $failure\n";
    }
}
my $duration = (time - $start) / 60;
print STDERR "$processed jobs processed in $duration minutes.  $goodFound were good.\n";
