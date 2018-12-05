=head1 Create a Matrix of Sample Contents

    build_sample_matrix.pl [options] outDir

This program takes the output of L<p3-list-reps.pl> when applied to the files created by L<get_sample_seeds.pl> and
produces a predictor directory (C<col.h>, C<row.h> C<X>) that characterizes each sample.

The standard input file to this method is tab-delimited with headers.  The first column contains an ID that consists of a sample name,
a period (C<.>), and a contig ID.  Only the sample name is used by this program.  The second column contains the representative
genome ID and the third column its name.  The final column is a similarity score.  A higher score indicates a better hit.

=head2 Parameters

The single positional parameter is the name of the output directory.  If the directory already exists, the C<--clear> option must
be specified.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

Additional command-line options are the following.

=over 4

=item clear

Erase the output directory before proceeding.  If this option is not specified, an error will occur if the output directory already
exists.

=item sim

The minimum permissible similarity score for a hit to be considered valid.  The default is C<100>.

=item nohead

If specified, the input file has no column headers.

=back

=cut

use strict;
use P3Utils;
use File::Copy::Recursive;

# Get the command-line options.
my $opt = P3Utils::script_opts('outDir', P3Utils::ih_options(),
        ['clear', 'erase output directory'],
        ['sim=i', 'minimum similarity score', { default => 100 }],
        ['nohead', 'input has no column headers']
        );
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (! -d $outDir) {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create $outDir: $!";
} elsif ($opt->clear) {
    # Erase the output directory.
    print "Erasing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not clear $outDir: $!";
} else {
    # Output directory already exists, but we are not clearing.
    die "$outDir already exists.";
}
# Get the similarity threshold.
my $minSim = $opt->sim;
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt, 1);
# This hash will track the genome IDs found.
my %genomes;
# This is a two-level hash keyed on {sampleID}{genomeID} that counts the number of occurrences.
my %hits;
# Loop through the input.
my ($count, $kept) = (0, 0);
while (! eof $ih) {
    my $line = <$ih>;
    my ($id, $rep, $name, $sim) = P3Utils::get_fields($line);
    # Only accept this line if the similarity is high enough.
    if ($sim >= $minSim) {
        # Remember the genome.
        if (! $genomes{$rep}) {
            print "Found $rep: $name.\n";
            $genomes{$rep} = $name;
        }
        # Count the hit.
        $kept++;
        my ($sampleID) = split /\./, $id;
        $hits{$sampleID}{$rep}++;
    }
    $count++;
    print "$count hits processed, $kept kept.\n" if $count % 1000 == 0;
}
# Now we have all the information we need.  Our first job is to create the col.h.  This assigns index numbers to the representative
# genomes.
my @cols = sort keys %genomes;
print "Writing column headers.";
open(my $oh, ">$outDir/col.h") || die "Could not open col.h: $!";
my $idx = 0;
for my $col (@cols) {
    print $oh join("\t", $idx, $col, $genomes{$col}) . "\n";
    $idx++;
}
print "  $idx genomes output.\n";
close $oh; undef $oh;
# Next we do the same with row.h, which assigns index numbers to the samples.
my @rows = sort keys %hits;
print "Writing row headers.";
open($oh, ">$outDir/row.h") || die "Could not open row.h: $!";
$idx = 0;
for my $row (@rows) {
    print $oh join("\t", $idx, $row) . "\n";
    $idx++;
}
print "  $idx samples output.\n";
close $oh; undef $oh;
# Finally, the matrix itself.
print "Writing matrix.";
open($oh, ">$outDir/X") || die "Could not open X: $!";
$kept = 0;
for my $r (@rows) {
    my @row = ();
    for my $c (@cols) {
        my $y = $hits{$r}{$c} // 0;
        push @row, "$y.0";
        $kept++ if $y;
    }
    print $oh join("\t", @row) . "\n";
}
print "  $kept cells were filled.\n";
close $oh; undef $oh;