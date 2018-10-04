=head1 Output Representative-Genome Groups

    p3-group-reps.pl [options] inDir

This script outputs the representative-genome groups present in a representative-genome directory.

The script operates on a representative server directory. The directory contains four files of paramount interest.

=over 4

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item rep_db.tbl

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID, (1) the ID of the genome's
representative, and (2) the similarity number.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

The output file will contain all of the genomes stored in the C<rep_db.tbl> file. The output file will be tab-delimited, with each record
containing (0) a representative genome ID, (1) a similarity score, (2) a represented genome ID, and (3) the represented genome's name.

=head2 Parameters

The positional parameter is the name of the input directory, which must contain the above four files.

The command-line options are as follows.

=over 4

=item sep

If specified, each representative-genome group will be separated by a line containing a double slash (C<//>).

=item verbose

If specified, progress messages will be written to the standard error output.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use Stats;

# Get the command-line options.
my $opt = P3Utils::script_opts('inDir',
        ['sep|group|g', 'separate representative-genome groups'],
        ['verbose|v', 'show progress on STDERR']
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Verify the parameters.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! -f "$inDir/6.1.1.20.fasta" || ! -f "$inDir/complete.genomes") {
    die "$inDir does not appear to contain a representative-genomes database.";
}
# Get the options.
my $verbose = $opt->verbose;
my $stats = Stats->new();
# Read in the representative-genome database.
my $repDB = RepGenomeDb->new_from_dir($inDir);
# Print the output header.
P3Utils::print_cols(['rep_id', 'score', 'genome_id', 'genome_name']);
# Loop through the representative genomes.
my $repList = $repDB->rep_list();
for my $repID (sort @$repList) {
    # Get this genome's rep object.
    my $repObject = $repDB->rep_object($repID);
    my $repName = $repObject->name;
    # Get the genome name and write out a line for this genome. The blank score only occurs for the representative itself.
    P3Utils::print_cols([$repID, '', $repID, $repName]);
    $stats->Add(groupOut => 1);
    # Get the list of represented genomes.
    my $genomesL = $repObject->rep_list();
    # Get the genome IDs sorted by similarity score.
    my @scoreList = sort { $b->[1] <=> $a->[1] } @$genomesL;
    my $memberCount = scalar @scoreList;
    print STDERR "Processing $repID: $repName. $memberCount members.\n" if $verbose;
    if ($memberCount > 0) {
        my @genomeIDs = map { $_->[0] } @scoreList;
        my $names = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name'], \@genomeIDs);
        my %nameHash = map { $_->[0] => $_->[1] } @$names;
        # Now output the results.
        for my $scoreTuple (@scoreList) {
            my ($genomeID, $score) = @$scoreTuple;
            my $name = $nameHash{$genomeID} // '<unknown>';
            P3Utils::print_cols([$repID, $score, $genomeID, $name]);
            $stats->Add(memberOut => 1);
        }
    } else {
        $stats->Add(nullGroup => 1);
    }
    # Optional separator prints here.
    if ($opt->sep) {
        print "//\n";
    }
}
print STDERR "All done.\n" . $stats->Show() if $verbose;
