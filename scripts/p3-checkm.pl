=head1 Run CheckM Against A List of PATRIC Genomes

    p3-checkm.pl [options] outDir

This script runs CheckM against a list of IDs for PATRIC genomes. The output directory will be used to contain working files.
In addition, the output lists will be put in there.

=head2 Parameters

The single positional parameter is the name of the output directory. This directory should not be shared with other
instances of this program, as it is used for working files.

The standard input can be overridden using the options in L<P3Utils/ih_options>. The input column should contain the genome
IDs.

Additional command-line options are the following.

=over 4

=item col

Index (1-based) of the column number to contain the key field. If a non-numeric value is specified, it is presumed
to be the value of the header in the desired column. The default is C<0>, which indicates the last column.

=item batchSize

Maximum number of lines to read in a batch. The default is C<10>.

=item nohead

Input file has no headers.

=back

=head2 Output Files

The genomes are presumed to have passed the prima facie check in L<p3-good-check.pl>. There are two main output files placed
in the output directory, both tab-delimited, each line consisting of a genome ID and name. The good genomes will be put in
C<good.patric.tbl> and the bad genomes in C<bad.patric.tbl>.

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use File::Copy::Recursive;
use SeedUtils;

# Get the command-line options.
my $opt = P3Utils::script_opts('outDir', P3Utils::ih_options(),
                ['col|c=s', 'column number (1-based) or name', { default => 0 }],
                ['batchSize|b=i', 'input batch size', { default => 10 }],
                ['nohead', 'file has no headers']
    );
# Create the statistics object.
my $stats = Stats->new();
# Insure we have an output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "No output directory specified.";
} elsif (-d $outDir) {
    print "Clearing $outDir.\n";
    File::Copy::Recursive::pathempty($outDir) || die "Could not clear output directory: $!";
} elsif (-f $outDir) {
    die "Invalid directory name $outDir.";
} else {
    print "Creating $outDir.\n";
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory: $!";
}
# Create the FASTA directory.
File::Copy::Recursive::pathmk("$outDir/fasta") || die "Could not create FASTA directory: $!";
# Create the work directory.
File::Copy::Recursive::pathmk("$outDir/checkm") || die "Could not create work directory: $!";
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Open the output files.
open(my $gh, "$outDir/good.patric.tbl") || die "Could not open good output file: $!";
open(my $bh, "$outDir/bad.patric.tbl") || die "Could not open bad output file: $!";
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # Get the names of these genomes.
    my $count = scalar @$couplets;
    $stats->Add(genomesIn => $count);
    print "Processing $count genomes in batch.\n";
    my $results = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name'], [map { $_->[0] } @$couplets]);
    # Loop through the genomes one at a time, creating the FASTA files and saving the names.
    my %names;
    for my $result (@$results) {
        my $id = $result->{genome_id};
        $names{$id} = $result->{genome_name};
        print "Reading sequences for $id: $names{$id}\n";
        $stats->Add(genomesProcessed => 1);
        my $triples = $p3->fasta_of($id);
        print "Creating FASTA file.\n";
        open(my $fh, ">$outDir/fasta/$id.fa") || die "Could not open output FASTA file: $!";
        for my $triple (@$triples) {
            print $fh ">$triple->[0] $triple->[1]\n$triple->[2]\n";
            $stats->Add(sequenceOut => 1);
        }
    }
    # Now we run CheckM on the FASTA files.
    my $cmd = "checkm lineage_wf --tmpdir $FIG_Config::temp -x fa --file $outDir/evaluate.txt $outDir/fasta $outDir/checkm";
    print "CheckM comment is: $cmd\n";
    SeedUtils::run($cmd);
    $stats->Add(runCheckM => 1);
    # Erase the FASTA files and CheckM work files.
    print "Clearing CheckM files.\n";
    File::Copy::Recursive::pathempty("$outDir/fasta") || die "Could not empty FASTA directory: $!";
    File::Copy::Recursive::pathempty("$outDir/checkm") || die "Could not empty checkM work directory: $!";
    # Read the evaluation file.
    print "Reading CheckM results.\n";
    open(my $ch, "<$outDir/evaluate.txt") || die "Could not open CheckM output file: $!";
    while (! eof $ch) {
        my $line = <$ih>;
        $stats->Add(checkMLineIn => 1);
        my @cols = split /\s+/, $line;
        if ($cols[1] =~ /(\d+\.\d+)/) {
            my $id = $1;
            $stats->Add(checkMResult => 1);
            my $score = $cols[13];
            my $contam = $cols[14];
            my $name = $names{$id} // "<unknown>";
            print "Score for $id ($name) is $score, contamination $contam.\n";
            if ($score >= 80 && $contam <= 10) {
                print "Genome is good.\n";
                $stats->Add(genomeGood => 1);
                print $gh "$id\t$name\n";
            } else {
                print "Genome is bad.\n";
                $stats->Add(genomeBad => 1);
                print $bh "$id\t$name\n";
            }
        }
    }
}
print "All done: \n" . $stats->Show();
