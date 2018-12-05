=head1 Extract Seed Proteins from Samples

    get_sample_seeds.pl [options] inDir outDir

This program extracts the seed proteins from processed binning samples and stores them in FASTA files in a specified
output directory.  The output files will be named C<samples.faa> and C<samples.fna> and will contain protein and DNA
sequences, respectively.

After running this script, use L<p3-list-reps.pl> to find the representative genomes for each protein, then
L<build_sample_matrix.pl> to format the output from it into sample characterizations.

=head2 Parameters

The positional parameters are the names of the input and output directories.

Additional command-line options are as follows.

=over 4

=item overwrite

If specified, then the output directory can already exist; otherwise, it must be created by this script.

=back

=cut

use strict;
use FIG_Config;
use P3Utils;
use Contigs;
use File::Copy::Recursive;


$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir',
        ['overwrite|old|clear', 'output directory can exist and will be overwritten']
        );
# Check the parameters.
my ($inDir, $outDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Input directory $inDir not found or invalid.";
} elsif (! $outDir) {
    die "No output directory specified.";
} elsif (-d $outDir && ! $opt->overwrite) {
    die "Output directory $outDir already exists.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
}
opendir(my $dh, $inDir) || die "Could not open bin directory: $!";
my @bins = grep { -s "$inDir/$_/bins.found.tbl" } readdir $dh;
closedir $dh;
my $found = scalar @bins;
print "$found bins found.\n";
my $count = 0;
open(my $oh, '>', "$outDir/samples.fna") || die "Could not open samples output: $!";
open(my $ph, '>', "$outDir/samples.faa") || die "Could not open samples output: $!";
for my $bin (@bins) {
    $count++;
    print "Reading contigs for $bin ($count of $found).\n";
    my $contigs = Contigs->new("$inDir/$bin/contigs.fasta", genomeID => $bin);
    open(my $ih, '<', "$inDir/$bin/bins.found.tbl") || die "Could not open $bin input: $!";
    my $subseq = 1;
    while (my $line = <$ih>) {
        chomp $line;
        my ($contig, $start, $dir, $len) = split /\t/, $line;
        if ($dir eq '-') {
            $start -= ($len - 1);
        }
        print "$subseq is $contig, $start $dir $len.\n";
        my $dna = $contigs->dna([$contig, $start, $dir, $len]);
        my $seq = $contigs->xlate([$contig, $start, $dir, $len]);
        print $ph ">$bin.$contig $subseq\n$seq\n";
        print $oh ">$bin.$contig $subseq\n$dna\n";
        $subseq++;
    }
    print "$subseq sequences output for $bin.\n";
}
close $oh;
