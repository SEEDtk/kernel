=head1 Check Seed Proteins in Downloaded Genomes

    p3x-check-downloaded-genomes.pl [options] gtoDir

This is a post-processing script for L<p3x-process-downloadable-genomes.pl>.  It presumes that the C<--gtos> option was used to
store L<GenomeTypeObject> files in the output directory.  The output file from the prior script will be read, and the putatively
good genomes examined.  The seed protein will be drawn out.  If it is of the proper size, it will be stored in an output FASTA file.
If it is not, the genome will be downgraded to not-good.  The output FASTA file may then be used as input to L<p3-check-reps.pl>.

=head2 Parameters

The single positional parameter is the name of the working directory containing the GTO files.  This will also be the output directory.
The updated input file will be C<genomes.tbl>.  The FASTA file will be C<seeds.fasta>.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  The input file should be the output from
L<p3x-process-downloadable-genomes.pl>.

=over 4

=item qual

If specified, then columns named C<completeness> and C<contamination> will be extracted, and used to compute the standard
quality measure (C<M> for contamination <= 5 and completeness >= 50, C<H> for contamination <= 5 and completeness >= 90, else C<L>).

=cut

use strict;
use P3DataAPI;
use P3Utils;
use GPUtils;
use Stats;
use GenomeTypeObject;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('gtoDir', P3Utils::ih_options(),
        ['qual', 'compute standard quality measure'],
        );
my $stats = Stats->new();
# Get the working directory.
my ($gtoDir) = @ARGV;
if (! $gtoDir) {
    die "No input directory specified.";
} elsif (! -d $gtoDir) {
    die "Input directory $gtoDir not found or invalid.";
}
# Get the options.
my $qual = $opt->qual;
print "Connecting to files.\n";
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read in the header line and find the critical columns.
my @qual;
if ($qual) {
    @qual = ('completeness', 'contamination');
}
my ($headers, $cols) = P3Utils::find_headers($ih, input => 'PATRIC ID', 'Good?', @qual);
my @outHeaders = @$headers;
if ($qual) {
    push @outHeaders, 'Quality';
}
# Open the output FASTA file and main output file.
open(my $fh, ">$gtoDir/seeds.fasta") || die "Could not open FASTA output: $!";
open(my $oh, ">$gtoDir/genomes.tbl") || die "Could not open main output: $!";
# Write the new headers.
P3Utils::print_cols(\@outHeaders, oh => $oh);
# Loop through the input.
my $count = 0;
while (! eof $ih) {
    # Get the next line.
    my $line = <$ih>;
    my @fields = P3Utils::get_fields($line);
    $stats->Add(lineIn => 1);
    $count++;
    # Get the genome ID and goodness flag.
    my ($genomeID, $goodFlag, $complete, $contam) = P3Utils::get_cols(\@fields, $cols);
    # If this genome is good, we need to process its seed protein.
    my $evalType;
    if (! $goodFlag) {
        $stats->Add(badGenome => 1);
        $evalType = 'bad';
    } else {
        $stats->Add(goodIn => 1);
        # Read the GTO and get its seed protein.
        print "Reading GTO for $genomeID ($count).\n";
        my $gto = GenomeTypeObject->create_from_file("$gtoDir/$genomeID.gto");
        my $seedProt = GPUtils::get_seed($gto);
        # Verify that the genome is genuinely good.
        if (! $seedProt) {
            $stats->Add(badSeed => 1);
            $fields[$cols->[1]] = '';
            $evalType = 'bad';
        } else {
            $stats->Add(goodGenome => 1);
            $evalType = 'good';
            # It's genuinely good, so write to the FASTA file.
            print $fh ">$genomeID\n$seedProt\n";
        }
    }
    # Compute the alternate quality measure.
    if ($qual) {
        my $qualType = 'L';
        if ($contam <= 5) {
            if ($complete >= 90) {
                $qualType = 'H';
            } elsif ($complete >= 50) {
                $qualType = 'M';
            }
        }
        $stats->Add("$evalType-$qualType" => 1);
        $stats->Add("genome-$qualType" => 1);
        push @fields, $qualType;
    }
    # Write the output record.
    P3Utils::print_cols(\@fields, oh => $oh);
    $stats->Add(lineOut => 1);
}
close $oh;
close $fh;
print "All done.\n" . $stats->Show();