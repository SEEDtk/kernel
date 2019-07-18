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
use EvalHelper;
use GenomeTypeObject;

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
# Loop through the gto directory.
print STDERR "Scanning $gtoDir.\n";
opendir(my $dh, $gtoDir) || die "Could not open $gtoDir: $!";
my @gtoFiles = grep { $_ =~ /\.gto$/ && -s "$gtoDir/$_" } readdir $dh;
closedir $dh;
my ($count, $total) = (0, scalar @gtoFiles);
print STDERR "$total GTOs found in $gtoDir.\n";
my %rolesBad;
for my $gtoFile (@gtoFiles) {
    # Read in the GTO.
    $count++;
    print STDERR "Processing $gtoFile: $count of $total.\n";
    my $gto = GenomeTypeObject->create_from_file("$gtoDir/$gtoFile");
    # Call the main processor.
    my $geo = EvalHelper::ProcessGto($gto, checkDir => $opt->checkdir, predictors => $opt->predictors,
        parallel => $opt->parallel, workDir => $opt->workdir, p3 => $p3);
    if ($opt->report) {
        my $roleReport = $geo->roleReport;
        for my $role (keys %$roleReport) {
            $rolesBad{$role}++;
        }
    }
    # Write the results.
    $gto->destroy_to_file("$gtoDir/$gtoFile");
}
if ($opt->report && $count > 0) {
    print STDERR "Generating role report.\n";
    my @roles = sort { $rolesBad{$b} <=> $rolesBad{$a} } keys %rolesBad;
    print "Role\tcount\tpercent\n";
    for my $role (@roles) {
        P3Utils::print_cols([$role, $rolesBad{$role}, $rolesBad{$role} * 100 / $count]);
    }
}