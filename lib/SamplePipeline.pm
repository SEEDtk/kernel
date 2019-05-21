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


package SamplePipeline;

    use strict;
    use warnings;
    use RASTlib;
    use Bin;
    use Bin::Analyze;
    use Bin::Improve;
    use Loader;
    use SeedUtils;
    use SeedAware;
    use RoleParse;

=head1 Process a Metagenome Sample Pipeline

This module will process a metagenome from the initial reads all the way through the expectation check. It must invoke
the following commands.

    spades.py -1 pair.1.fastq -2 pair.2.fastq -s singles.fastq -o Assembly --meta
    bins_coverage Assembly/contigs.fasta .
    bins_generate .
    ... call RastBins ...
    p3x-eval-bins .
    bins_expect_check --genus --input=expect.tbl .

This process assembles the reads into contigs, separates them into bins containing the major populations, sorted by genus,
submits them to RAST to compute the proteins, and then analyzes the results.

The four input files are

=over 4

=item pair.1.fastq

The fastq file for the first set of paired reads. For an HMP microbiome, this has the name
I<sample>C<.denovo_duplicates_marked.trimmed.1.fastq>.

=item pair.2.fastq

The fastq file for the second set of paired reads. For an HMP microbiome, this has the name
I<sample>C<.denovo_duplicates_marked.trimmed.2.fastq>.

=item singles.fastq

The fastq file for the singleton reads. For an HMP microbiome, this has the name
I<sample>C<.denovo_duplicates_marked.trimmed.singleton.fastq>.

=item expect.tbl

A file listing the species expected in the metagenome. The first column should be the
species name and the third column the coverage depth. For an HMP microbiome, this has the
name I<sample>C<_abundance_table.tsv>.

=back

These four files are specified as parameters to this object's main method. If one of the FASTQ
files is missing, it is not provided to the spades assembly script. If the expectation table is
missing, the expectation report is not run.

The directory in which the program runs will be filled with many different files. The following are
of interest to the Shrub database.

=over 4

=item binX.gto

A L<GenomeTypedObject>, stored in JSON format, for one of the computed bins, where I<X> is the bin number.

=item expect.report.txt

The expectation report, describing the genuses expected and found in the sample.

=item bins.rast.json

The bins themselves, in JSON format, with the universal proteins computed.

=item bins.report.txt

The bins report, containing statistics about the bins and their contents.

=item ref.genome.score.tbl

The primer protein BLAST hits in the contigs computed by the assembly.

=back

=head2 Public Methods

=head3 Process

    Process($workDir, %options);

Process and analyze a metagnome sample. The sample will be assembled into contigs and separated into bins, then
the bins will be processed by RAST. Optionally, the results can be compared to expected results.

=over 4

=item workDir

The working directory into which all output files should be placed.

=item options

A hash containing the following options.

=over 4

=item f1

The name of the fastq file for a first set of paired reads.

=item f2

The name of the fastq file for a second set of paired reads.

=item fs

The name of the fastq file for the singleton reads.

=item expect

The name of a tab-delimited file containing the expected species population.

=item eNameCol

The index (1-based) of the expectation file column containing the species names. The default is C<1>.

=item eDepthCol

The index (1-based) of the expectation file column containing the species coverage depth. The default is C<3>.

=item force

If TRUE, then all the files will be regenerated, regardless of whether or not they already exist. If FALSE, the existence
of files will presuppose that a run was interrupted and needs to be finished.

=item large

If TRUE, the memory footprint will be made larger, but the assembly process will be slower. Use this for larger read sets.

=item user

The apprioriate RAST user name. The default is taken from the C<RASTUSER> environment variable.

=item password

The appropriate RAST password. The default is taken from the C<RASTPASS> environment variable.

=item sleep

The number of seconds to sleep between polling requests to the RAST server. The default is C<60>.

=item engine

Type of binning engine to use-- C<s> for standard or C<2> for alternate.

=item noIndex

If TRUE, the annotated genomes produced will not be indexed in PATRIC.

=item uniRoles

The name of a file containing the universal roles.  Each record of the file should contain a role ID, a checksum,
and a role description.  The default is C<uniRoles.tbl> in the SEEDtk global directory.

=back

=back

If none of the three FASTQ names are specified, then no assembly step will be performed.

The assembled contigs must exist in C<contigs.fasta> in the working directory and the coverage information
in the same directory's C<output.contigs2reads.txt> file. If no expectation file is given, no expectation
report will be produced.

=cut

sub Process {
    my ($workDir, %options) = @_;
    # Extract the options.
    my $f1q = $options{f1};
    my $f2q = $options{f2};
    my $fsq = $options{fs};
    my $engine = $options{engine} // 's';
    my $eNameCol = $options{eNameCol} // 1;
    my $eDepthCol = $options{eDepthCol} // 3;
    my %rastOptions = (user => $options{user}, password => $options{password}, 'sleep' => $options{sleep}, noIndex => $options{noIndex});
    my $force = $options{force} // 0;
    my $uFile = $options{uniRoles} // "$FIG_Config::global/uniRoles.tbl";
    # Build the uni-role hash.
    my %uniRoles;
    open(my $uh, '<', $uFile) || die "Could not open universal role file: $!";
    while (! eof $uh) {
        my $line = <$uh>;
        chomp $line;
        my ($role, $check, $name) = split /\t/, $line;
        $uniRoles{$check} = [$role, $name];
    }
    close $uh;
    print scalar(%uniRoles) . " universal roles found.\n";
    # Do we need to assemble the contigs?
    my $assemble = 0;
    my $haveContigs = (-s "$workDir/contigs.fasta" && -s "$workDir/output.contigs2reads.txt");
    print "Contigs " . ($haveContigs ? '' : 'not ') . "found in directory.\n";
    if ($f1q || $f2q || $fsq) {
        if ($force || ! $haveContigs) {
            $assemble = 1;
            my $fileCount = 0;
            for my $file ($f1q, $f2q, $fsq) {
                if ($file && -f $file) {
                    $fileCount++;
                }
            }
            if (! $fileCount) {
                die "Contigs not assembled, but no read files found.";
            }
        }
    } elsif (! $haveContigs) {
        die "Contigs not assembled, and no fastq files specified.";
    }
    if ($assemble) {
        # We need to assemble.
        my $cmd = "spades.py";
        my ($threads, $mem) = (8, 250);
        if ($options{large}) {
            ($threads, $mem) = (4, 500);
        }
        # Add the directory and the style parms.
        my @parms = ('-o', "$workDir/Assembly", '--meta', '--threads', $threads, '--memory', $mem);
        # Determine the type of reads. If there is only "s", we do interleaved.
        if (! $f1q && ! $f2q) {
            push @parms, '-12', $fsq;
        } else {
            push @parms, '-1', $f1q, '-2', $f2q;
            if ($fsq) {
                push @parms, '-s', $fsq;
            }
        }
        # Find the command.
        my $cmdPath = SeedAware::executable_for($cmd);
        die "Could not find $cmd." if ! $cmdPath;
        # Execute the command.
        print "Running spades: " . join(" ", @parms) . "\n";
        my $rc = system($cmdPath, @parms);
        die "Error exit $rc from $cmd." if $rc;
        # Now run the coverage computation.
        $rc = system('bins_coverage', "$workDir/Assembly/contigs.fasta", $workDir);
        die "Error exit $rc from bins_coverage." if $rc;
        $force = 1;
    }
    # At this point, we have the contigs and the coverage data.
    if ($force || ! -s "$workDir/bins.json") {
        # We need to generate bins.
        my $command = join('_', "bin$engine", 'generate');
        print "Generating bins with $command.\n";
        my $rc = system($command, $workDir);
        die "Error exit $rc from bins_generate." if $rc;
        $force = 1;
    }
    # At this point, we have the bins JSON file.
    if ($force || ! -f "$workDir/bins.report.txt") {
        # We need to run RAST on the bins. This is done internally, and we need some parameters.
        # First, we set the "partial" option to the inverse of "force".
        $rastOptions{partial} = ($force ? 0 : 1);
        # Get the files.
        my $binJsonFile = "$workDir/bins.json";
        my $contigFastaFile = "$workDir/sample.fasta";
        RastBins(\%uniRoles, $binJsonFile, $contigFastaFile, $workDir, %rastOptions);
        $force = 1;
    }
    # Now we need to insure we have an evaluation.
    if ($force || ! -s "$workDir/Eval/index.tbl") {
        if (! -s "$workDir/bins.rast.json") {
            # No bins to evaluate.
            mkdir "$workDir/Eval" || die "Could not create Eval for $workDir: $!";
            open(my $oh, '>', "$workDir/Eval/index.tbl") || die "Could not open eval output file: $!";
            print $oh "No bins found.\n";
        } else {
            # Configure the options.  We always do a deep analysis, but we need to know if the genome is indexed.
            my @options = ('--deep');
            if ($options{noIndex}) {
                push @options, '--noIndex';
            }
            my $rc = system('p3x-eval-bins', @options, $workDir);
            die "Error exit $rc from p3x-eval-bins." if $rc;
        }
    }
}

=head3 RastBins

    SamplePipeline::RastBins($shrub, $binJsonFile, $contigFastaFile, $workDir, %options);

Run the bins through the PATRIC RAST to compute the universal roles and create genome estimates. This method
will output a json-format L<GenomeTypedObject> for each bin, with the name C<binX.gto>, with I<X> being the
bin nunber. It will create a new JSON file for the bins themselves-- C<bins.rast.json>-- with universal role
information included. Finally, it will produce an analysis report in C<bins.report.txt> describing the quality
of the bins.

=over 4

=item uniRoleH

A hash mapping each universal role's checksum to a 2-tuple of (0) ID and (1) description.

=item binJsonFile

The name of the file containing the bin objects in JSON format.

=item contigFastaFile

The name of a FASTA file containing the assembled contigs being binned.

=item workDir

The name of the working directory to contain the output files.

=item options

A hash containing zero or more of the following options.

=over 8

=item user

The appropriate RAST user name. The default is taken from the C<RASTUSER> environment variable.

=item password

The appropriate RAST password. The default is taken from the C<RASTPASS> environment variable.

=item sleep

The number of seconds to sleep between polling requests to the RAST server. The default is C<60>.

=item partial

If TRUE, only bins that have not yet been submitted to RAST will be processed. Otherwise, all bins will
be processed.

=item RETURN

Returns a reference to a list of L<Bin> objects representing the bins found.

=back

=back

=cut

sub RastBins {
    my ($uniRoleH, $binJsonFile, $contigFastaFile, $workDir, %options) = @_;
    # Get the options.
    my $rastUser = $options{user} // $ENV{RASTUSER};
    my $rastPass = $options{password} // $ENV{RASTPASS};
    my $sleep = $options{'sleep'} // 60;
    my $partial = $options{partial} // 0;
    # Count the universal roles.
    my $totUnis = scalar keys %$uniRoleH;
    # Build the uni-role hashes.
    my (%uniCheck, %uniName);
    for my $check (keys %$uniRoleH) {
        my $tuple = $uniRoleH->{$check};
        $uniCheck{$check} = $tuple->[0];
        $uniName{$tuple->[0]} = $tuple->[1];
    }
    # Get the loader object.
    my $loader = Loader->new();
    my $stats = $loader->stats;
    # Create the bin improvement object.
    my $improver = Bin::Improve->new($workDir, stats => $stats);
    # Initialize the role hashes.
    print "Reading role files.\n";
    my ($nMap, $cMap) = EvalCon::LoadRoleHashes("$FIG_Config::p3data/roles.in.subsystems", $stats);
    # Create the GEO options.
    my $p3 = P3DataAPI->new();
    my %gOptions = (roleHashes => [$nMap, $cMap], detail => 2, p3 => $p3, stats => $stats);
    # Read in the bins.
    print "Reading bins from $binJsonFile.\n";
    my $binList = Bin::ReadBins($binJsonFile);
    # Create the RAST option hash.
    my %rastOpts = (user => $rastUser, password => $rastPass, 'sleep' => $sleep);
    # Loop through the bins, processing them one at a time. For each bin, we read the whole sample
    # file to get the contigs. Only one bin's worth of contigs is kept in memory at a time, at the cost of
    # a lot of extra reading.
    my $binNum = 0;
    for my $bin (@$binList) {
        my $taxonID = $bin->taxonID;
        my $name = $bin->name;
        $stats->Add(bins => 1);
        $binNum++;
        print "Processing bin $binNum - $taxonID: $name.\n";
        my $binFile = "$workDir/bin$binNum.gto";
        my $gto;
        if ($partial && -s $binFile) {
            print "Reading bin GTO from $binFile.\n";
            $stats->Add(binsSkipped => 1);
            $gto = SeedUtils::read_encoded_object($binFile);
        } else {
            my %contigs = map { $_ => 1 } $bin->contigs;
            # Now we read the sample file and keep the contig triples.
            my $triples = [];
            my $ih = $loader->OpenFasta(sampleContigs => $contigFastaFile);
            my $binFastaFile = "$workDir/bin$binNum.fa";
            open(my $oh, '>', $binFastaFile) || die "Could not open FASTA output file for bin $binNum.";
            my $triple = $loader->GetLine(sampleContigs => $ih);
            while (defined $triple) {
                my $contigID = $triple->[0];
                if ($contigs{$contigID}) {
                    $stats->Add(contigsKept => 1);
                    push @$triples, $triple;
                    print $oh ">$triple->[0] $triple->[1]\n$triple->[2]\n";
                }
                $triple = $loader->GetLine(sampleContigs => $ih);
            }
            close $oh;
            my $contigCount = scalar @$triples;
            print "Submitting $binNum to RAST: $contigCount contigs.\n";
            $gto = RASTlib::Annotate($triples, $taxonID, $name, %rastOpts);
            # Check to see if we can improve this bin.
            if ($improver->eligible($gto)) {
                # Yes.  Try to improve it.
                print "Attempting to improve $binNum.\n";
                $triples = $improver->Process($bin, $gto, $binFastaFile, $triples);
                if ($triples) {
                    print "Submitting improved bin to RAST.\n";
                    $gto = RASTlib::Annotate($triples, $taxonID, "$name cleaned", %rastOpts);
                }
            }
            print "Spooling genome to $workDir.\n";
            SeedUtils::write_encoded_object($gto, "$workDir/bin$binNum.gto");
        }
        print "Fixing bin object.\n";
        # Clear the bin's current universal protein list.
        $bin->replace_prots();
        # Search the genome for universal roles.
        my $flist = $gto->{features};
        for my $feature (@$flist) {
            my $fun = $feature->{function};
            if ($fun) {
                my @roles = SeedUtils::roles_of_function($fun);
                for my $role (@roles) {
                    my $checksum = RoleParse::Checksum($role);
                    my $roleID = $uniCheck{$checksum};
                    if ($roleID) {
                        $bin->incr_prot($roleID, 1);
                    }
                }
            }
        }
        # Adjust the contig list.
        my @cList = map { $_->{id} } @{$gto->{contigs}};
        $bin->AdjustContigList(\@cList);
        print "Final bin size is " . scalar(@cList) . " contigs with length " . $bin->len . ".\n";
    }
    # Output the new bins. Note we don't sort any more because we need to preserve the
    # bin number.
    my $binOutputFile = "$workDir/bins.rast.json";
    print "Spooling bins to $binOutputFile.\n";
    open(my $oh, ">$binOutputFile") || die "Could not open bins.rast.json output file: $!";
    for my $bin (@$binList) {
        $bin->Write($oh);
    }
    close $oh;
    # Create an analyzer object.
    my $analyzer = Bin::Analyze->new(totUnis => $totUnis, minUnis => (0.8 * $totUnis));
    # Read the reference genome file. We need this for the report.
    my $refScoreFile = "$workDir/ref.genomes.scores.tbl";
    if (-s $refScoreFile) {
        my %genomes;
        print "Reading reference genomes from $refScoreFile.\n";
        my $ih = $loader->OpenFile(refGenomes => $refScoreFile);
        while (my $refFields = $loader->GetLine(refGenomes => $ih)) {
            my ($contig, $genome, $score, $name) = @$refFields;
            $genomes{$genome} = $name;
        }
        $analyzer->SetGenomes(\%genomes);
    }
    # Write the report.
    open(my $rh, ">$workDir/bins.report.txt") || die "Could not open report file: $!";
    $analyzer->BinReport($rh, \%uniName, \@$binList);
    close $rh;
}


1;