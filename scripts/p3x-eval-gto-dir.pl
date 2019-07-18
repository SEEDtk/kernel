=head1 Evaluate All GenomeTypeObjects in a Directory

    p3x-eval-gto-dir.pl [options] gtoDir

This script evaluates all L<GenomeTypeObject> files in a directory and stores the quality data back into them.
Progress messages will be written to STDERR.

=head2 Parameters

The positional parameter is the name of the L<GenomeTypeObject> directory.

Additional command-line options are as follows:

=over 4

=item checkDir

The name of the directory containing the reference genome table and the completeness data files. The default
is C<CheckG> in the SEEDtk global data directory.

=item predictors

The name of the directory containing the role definition files and the function predictors for the consistency
checking. The default is C<FunctionPredictors> in the SEEDtk global data directory.

=item parallel

The number of parallel processes to run when applying the function predictors. The default is C<8>.

=item workDir

Name of a working directory for creating the matrix.  If none is specified, a temporary directory will be used.

=item report

If specified, a report of the most common problematic roles will be written to the standard output.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use EvalCon;
use EvalCom::Rep;
use EvalCom::Tax;
use GenomeTypeObject;
use GEO;

# Get the command-line options.
my $opt = P3Utils::script_opts('gtoFile outFile outHtml',
        ['ref|r=s', 'reference genome ID (implies deep)'],
        ['deep', 'if specified, the genome is compared to a reference genome for more detailed analysis'],
        ['checkDir=s', 'completeness data directory', { default => "$FIG_Config::p3data/CheckG" }],
        ['predictors=s', 'function predictors directory', { default => "$FIG_Config::p3data/FunctionPredictors" }],
        ['parallel=i', 'parallelism to use in matrix evaluation', { default => 8 }],
        ['workDir=s', 'name of a working directory for the evaluation matrix'],
        ['report', 'produce a problematic roles report on the standard output']
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the input parameters.
my ($gtoDir) = @ARGV;
if (! $gtoDir) {
    die "No input GTO file specified.";
} elsif (! -d $gtoDir) {
    die "Invalid or missing GTO directory $gtoDir.";
}
# Create the work directory.
my $workDir = $opt->workdir;
my $tmpObject;
if (! $workDir) {
    $tmpObject = File::Temp->newdir();
    $workDir = $tmpObject->dirname;
} elsif (! -d $workDir) {
    File::Copy::Recursive::pathmk($workDir) || die "Could not create work directory: $!";
}
# Create the consistency helper.
my $evalCon = EvalCon->new(predictors => $opt->predictors);
# Get access to the statistics object.
my $stats = $evalCon->stats;
# Create the completeness helper.
my $checkDir = $opt->checkdir // "$FIG_Config::p3data/CheckG";
my ($nMap, $cMap) = $evalCon->roleHashes;
my %evalOptions = (logH => \*STDERR, stats => $stats);
my $evalG;
if (-s "$checkDir/REP") {
    open(my $xh, '<', "$checkDir/REP") || die "Could not open REP file: $!";
    my $k = <$xh>;
    chomp $k;
    $evalG = EvalCom::Rep->new($checkDir, %evalOptions, K => $k);
} else {
    $evalG = EvalCom::Tax->new($checkDir, %evalOptions, roleHashes=> [$nMap, $cMap]);
}
my %geoOptions = (roleHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, detail => 1);
# Loop through the gto directory.
print STDERR "Scanning $gtoDir.\n";
opendir(my $dh, $gtoDir) || die "Could not open $gtoDir: $!";
my @gtoFiles = grep { $_ =~ /\.gto$/ && -s "$gtoDir/$_" } readdir $dh;
closedir $dh;
my ($count, $total) = (0, scalar @gtoFiles);
print STDERR "$total GTOs found in $gtoDir.\n";
my %rolesBad;
my @gtoList;
for my $gtoFile (@gtoFiles) {
    # Collect the GTO files.
    $count++;
    print STDERR "Loading $gtoFile: $count of $total.\n";
    push @gtoList, $gtoFile;
    if (scalar @gtoList >= 50) {
        ProcessGtoList(\@gtoList, \%rolesBad);
        @gtoList = ();
    }
}
if (scalar @gtoList) {
    ProcessGtoList(\@gtoList, \%rolesBad)
}
if ($opt->report && $count > 0) {
    print STDERR "Generating role report.\n";
    my @roles = sort { $rolesBad{$b} <=> $rolesBad{$a} } keys %rolesBad;
    print "Role\tcount\tpercent\n";
    for my $role (@roles) {
        P3Utils::print_cols([$role, $rolesBad{$role}, $rolesBad{$role} * 100 / $count]);
    }
}

sub ProcessGtoList {
    my ($gtoList, $rolesBad) = @_;
    # Create the eval matrix for the consistency checker.
    $evalCon->OpenMatrix($workDir);
    # Loop through the GTOs, gathering completeness data and preparing for the consistency check.
    my %geoMap;
    my %gtoMap;
    for my $gtoFile (@$gtoList) {
        my $gto = GenomeTypeObject->create_from_file("$gtoDir/$gtoFile");
        my $geo = GEO->CreateFromGto($gto, %geoOptions);
        my $genomeID = $gto->{id};
        $gtoMap{$gtoFile} = $gto;
        $geoMap{$genomeID} = $geo;
        # Open the output file for the quality data.
        my $qFile = "$workDir/$genomeID.out";
        open(my $oh, '>', $qFile) || die "Could not open work file: $!";
        # Output the completeness data.
        $evalG->Check2($geo, $oh);
        close $oh;
        $evalCon->AddGeoToMatrix($geo);
    }
    $evalCon->CloseMatrix();
    # Evaluate the consistency.
    my $rc = system('eval_matrix', "-p", $opt->parallel, '-q', $evalCon->predictors, $workDir, $workDir);
    if ($rc) {
        die "EvalCon returned error code $rc.";
    }
    # Store the quality metrics in the GEOs.
    for my $genomeID (keys %geoMap) {
        my $qFile = "$workDir/$genomeID.out";
        $geoMap{$genomeID}->AddQuality($qFile);
    }
    # Update the GTOs from the GEOs.
    for my $gtoFile (@$gtoList) {
        my $gto = $gtoMap{$gtoFile};
        my $genomeID = $gto->{id};
        my $geo = $geoMap{$genomeID};
        $geo->UpdateGTO($gto);
        $gto->destroy_to_file("$gtoDir/$gtoFile");
        # Update the report if necessary.
        if ($opt->report) {
            my $roleReport = $geo->roleReport;
            for my $role (keys %$roleReport) {
                $rolesBad->{$role}++;
            }
        }
    }
}