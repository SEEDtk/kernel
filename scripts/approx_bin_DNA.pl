use SeedUtils;
use strict;
use Data::Dumper;
use RepKmers;
use FastA;
use FastQ;
use ScriptUtils;
use Stats;
use Signatures;
use File::Copy::Recursive;

=head1 Use Genome Signatures to Bin DNA

    approx_bin_DNA outDir file1 file2 ... fileN <signatures

Apply signatures from L<aa_signatures_of_genomes.pl> to reads or contigs in order to form
bins. Each input sequence will be binned according to which genome's signatures are most
common.

=head2 Parameters

The standard input is a two-column file of protein signatures, tab-delimited, with the first
column being the signatures and the second column the target representative genome ID.

The positional parameters are an output directory name followed by the names of the input files.
If the input is paired-end reads, the paired files must be placed together in the sequence
I<left1 right1 left2 right2 ... leftN rightN>.

The command-line options are those found in L<ScriptUtils::ih_options> (to specify the standard
input), L<Shrub::script_options> to specify the database, and the following.

=over 4

=item fasta

If specified, the sequence input is presumed to be FASTA files.

=item interlaced

If specified, the sequence input is presumed to be interlaced reads.

=item paired

If specified, the sequence input is presumed to be paired-end reads. This is the default.

=item min

Minimum signature score necessary for a sequence to be binned. A sequence is put in the bin for which it
has the highest score, but only if that score is greater than this number. The default is C<3>.

=item scoring

Type of scoring algorithm. Currently, the permissible values are

=over 8

=item raw

Score by the raw number of signature occurrences times 1000 divided by the number of kmers in the sequence.

=item scaled

Score by the number of signature occurrences times 1000 divided by the number kmers in the sequence times the total
number of signatures for the genome.

=back

=item geneticCode

Genetic code to use for translating DNA into proteins. The default is C<11>.

=back

=head2 Ouput

Executing the command produces a directory (the data directory specified in the first positional parameter) in
which there is a file for each genome with signatures.

Each output file will have the same name as the representative genome ID identified by the dominant signatures,
with a suffix of C<.fasta> or C<.fastq> depending on the input format. The fastq file will be interlaced.

=cut

use constant CLASSES => { raw => 'Signatures::Raw', scaled => 'Signatures::Scaled' };
$| = 1;
my $opt = ScriptUtils::Opts('outDir file1 file2 ... fileN',
        ScriptUtils::ih_options(),
        ['mode', 'hidden', { one_of => [['paired', 'paired-end FASTQ input'],
                                        ['fasta', 'FASTA sequence input'],
                                        ['interlaced', 'interlaced FASTQ input']],
                             default => 'paired' }],
        ['min=f', 'minimum score for binning', {default => 3}],
        ['scoring=s', 'scoring algorithm', { default => 'raw' }],
        ['geneticCode=i', 'protein translation code', { default => 11 }],
        );
my ($outDir, @files) = @ARGV;
# Verify the output directory.
if (! -$outDir || -f $outDir) {
    die "No output directory specified or output directory is a file.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
}
# Compute the signature class.
my $class = CLASSES->{lc $opt->scoring};
die "Invalid scoring option " . $opt->scoring . "." if (! $class);
# Include the scoring module.
my $requireFile = join('/', split /::/, $class). ".pm";
require $requireFile;
# Create the statistics object.
my $stats = Stats->new();
# Get the input signatures.
print "Reading input signatures.\n";
my $ih = ScriptUtils::IH($opt->input);
my $sigsObject = Signatures::new($class, $ih, $opt->min, $stats, geneticCode => $opt->geneticcode);
# This hash will contain the output file handles.
my %outputFiles;
# Compute the input mode.
my $mode = $opt->mode;
# Compute the file name suffix.
my $suffix = ($mode eq 'fasta' ? 'fasta' : 'fastq');
# This will contain the current reader object.
my $reader;
# Compute the input file list. If we are using paired files, this pairs them together.
my $filesList = FastQ::OrganizeFiles(($mode ne 'paired'), @files);
my $fileTotal = scalar @$filesList;
my $fileCount = 0;
# Loop through the input file sets.
for my $fileSet (@$filesList) {
    $fileCount++;
    print "Processing file set $fileCount of $fileTotal.\n";
    if ($mode eq 'fasta') {
        $reader = FastA->new(@$fileSet);
    } else {
        $reader = FastQ->new(@$fileSet);
    }
    # Loop through the sequences.
    while ($reader->next) {
        $stats->Add(inputRecords => 1);
        my @bins = $sigsObject->ComputeBin($reader);
        # Output the sequence to the bins.
        for my $bin (@bins) {
            my $fh = $outputFiles{$bin};
            if (! $fh) {
                open($fh, ">$outDir/$bin.$suffix") || die "Could not open output file for $bin: $!";
                $stats->Add(outputBins => 1);
                $outputFiles{$bin} = $fh;
            }
            $reader->Write($fh);
        }
    }
    $stats->Add(inputSets => 1);
}
print "All done.\n" . $stats->Show();