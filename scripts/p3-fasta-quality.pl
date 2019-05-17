=head1 Annotate and Quality-Check a FASTA File

    p3-fasta-quality.pl [options] taxonID name outDir

This script will annotate a FASTA file and perform quality analysis. The user must be logged into RAST using L<p3-login.pl>.
The output will be a GTO named C<bin.gto> in the specified output directory containing the annotated genome. There will
be subdirectories C<EvalG> and C<EvalCon> containing the output files from the two quality-check tools. The file
C<report.html> will contain a full quality report.

=head2 Parameters

The positional parameters are a taxonomic or genome ID in which the genome is believed to reside, the name to give the
genome, and the name of the output directory.

The standard input should contain the FASTA file and can be overridden using the options in L<P3Utils/ih_options>.

Additional command-line options are the following.

=over 4

=item domain

The domain of the new genome-- C<B> for bacteria, C<A> for archaea, and so forth. The default is
C<B>.

=item geneticCode

The genetic code of the new genome. The default is C<11>.

=item predictors

Name of the function predictors directory. The default is FunctionPredictors in the SEEDtk global data directory.
If this option is specified, the role files in the predictors directory will be used instead of the global role files.

=item sleep

Sleep interval in seconds while waiting for the RAST job to complete. The default is C<60>.

=item tDir

The name of the directory containing the web page templates.

=back

=cut

use strict;
use P3Utils;
use RASTlib;
use BinningReports;
use gjoseqlib;
use File::Copy::Recursive;
use GenomeTypeObject;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('taxonID name outDir', P3Utils::ih_options(),
        ["predictors=s", 'function predictors directory'],
        ["domain|d=s", "domain (A or B) of the new genome", { default => 'B' }],
        ["geneticCode=i", "genetic code for the new genome", { default => 11 }],
        ["sleep=i", "sleep interval for status polling", { default => 60 }],
        ['tDir|tdir|templates=s', 'template file directory', { default => "$FIG_Config::mod_base/kernel/lib/BinningReports" }],
        );
# Get the positional parameters.
my ($taxonID, $name, $outDir) = @ARGV;
if (! $taxonID) {
    die "No genome or taxon ID specified.";
} elsif (! $name) {
    die "No genome name specified.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory: $!";
} else {
    print "Clearing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not erase output directory: $!";
}
print "Processing options.\n";
my ($rolesToUse, $roleFile, $predictors);
if ($opt->predictors) {
    $predictors = $opt->predictors;
    $rolesToUse = "$predictors/roles.to.use";
    $roleFile = "$predictors/roles.in.subsystems";
} else {
    $predictors = "$FIG_Config::p3data/FunctionPredictors";
    $rolesToUse = "$FIG_Config::p3data/roles.to.use";
    $roleFile = "$FIG_Config::p3data/roles.in.subsystems";
}
print "EvalCon files are $predictors, $rolesToUse, and $roleFile.\n";
my $domain = $opt->domain;
my $geneticCode = $opt->geneticcode;
# Get the contigs from the file. We form the contigs into a FASTA string.
print "Reading FASTA data.\n";
# Open the input file.
my $ih = P3Utils::ih($opt);
my $contigs = gjoseqlib::read_fasta($ih);
print scalar(@$contigs) . " contigs found.\n";
print "Annotating genome.\n";
# Invoke the RAST service.
my $gto = RASTlib::Annotate($contigs, $taxonID, $name, user => undef, password => undef,
        domain => $domain, geneticCode => $geneticCode, sleep => $opt->sleep);
my $genomeID = $gto->{id};
print "Writing $genomeID to $outDir/bin.gto.\n";
SeedUtils::write_encoded_object($gto, "$outDir/bin.gto");
# Run the quality checks.
my $cmd = "check_gto --eval --quiet $outDir/EvalG $outDir/bin.gto";
print "Running $cmd.\n";
SeedUtils::run($cmd);
$cmd = "gto_consistency $outDir/bin.gto $outDir/EvalCon $predictors $roleFile $rolesToUse";
print "Running $cmd.\n";
SeedUtils::run($cmd);
# Create the role maps.
print "Reading role file.\n";
my (%cMap, %nameMap);
open(my $rh, '<', $roleFile) || die "Could not open role file: $!";
while (! eof $rh) {
    my $line = <$rh>;
    my ($id, $cksum, $name) = P3Utils::get_fields($line);
    $cMap{$cksum} = $id;
    $nameMap{$id} = $name;
}
# Prepare for the quality report.
print "Indexing GTO.\n";
GenomeTypeObject->initialize($gto);
print "Producing quality report.\n";
BinningReports::UpdateGTO($gto, "$outDir/EvalCon", "$outDir/EvalG", \%cMap);
$gto->destroy_to_file("$outDir/bin.gto");
print "Producing HTML.\n";
# Load the detail template.
my $tDir = $opt->tdir;
open(my $th, "<$tDir/details.tt") || die "Could not open detail template file: $!";
my $detailsT = join("", <$th>);
close $th; undef $th;
# Build the HTML prefix and suffix.
my $prefix = "<html><head>\n<style type=\"text/css\">\n";
open($th, "<$tDir/packages.css") || die "Could not open style file: $!";
while (! eof $th) {
    $prefix .= <$th>;
}
close $th; undef $th;
$prefix .= "</style></head><body>\n";
my $suffix = "\n</body></html>\n";
my $html = BinningReports::Detail({}, undef, \$detailsT, $gto, \%nameMap);
open(my $oh, ">$outDir/report.html");
print $oh "$prefix\n$html\n$suffix\n";
close $oh;
