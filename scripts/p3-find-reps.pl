=head1 Find Representative Genomes

    p3-find-reps.pl [options] inDir

This script looks at an input list of PATRIC and/or Shrub genome IDs and outputs the representative genome for each one.

The script operates on a representative server directory. The directory contains three files of paramount interest.

=over 4

=item 6.1.1.20.fasta

A FASTA file containing the identifying protein (Phenylalanyl tRNA synthetase alpha chain) for each representative genome. The
genome ID should be the sequence ID; however, if a comment is present, it will be assumed the sequence ID is a feature ID and
the comment contains the genome ID.

=item complete.genomes

A tab-delimited file (no headers) containing one record per genome with the columns (0) the genome ID and (1) the genome name.

=item K

A parameter file containing two records. The first record contains the protein kmer size (e.g. C<8>) and the second contains the minimum
similarity number for a genome to be considered represented (e.g. C<100>).

=back

=head2 Parameters

The positional parameters are the input directory and the output directory. The input directory must contain the above three files.

The standard input can be overridden using the options in L<P3Utils/ih_options>. It should contain PATRIC genome IDs in the key column.

Additional command-line options are those given in L<P3Utils/col_options> (to select the input column) and L<Shrub/script_options>
(to identify the L<Shrub> database) plus the following.

=over 4

=item fasta

If specified, the standard input will be a FASTA file of protein sequences with genome or feature IDs.  If this option is
specified, the databases are not accessed; instead, the protein sequences are read from the input.  The FASTA comment will
be used in place of the genome name.

=item verbose

If specified, progress messages will be written to STDERR.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use File::Copy::Recursive;
use Stats;
use RoleParse;
use FastA;
use Shrub;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('inDir outDir', P3Utils::col_options(), P3Utils::ih_options(), Shrub::script_options(),
        ['fasta', 'input is FASTA proteins, not genome IDs'],
        ['verbose|debug|v', 'write progress messages to STDERR']
        );
my $debug = $opt->verbose;
# Create the statistics object.
my $stats = Stats->new();
# Verify the parameters.
my ($inDir) = @ARGV;
if (! $inDir) {
    die "No input directory specified.";
} elsif (! -d $inDir) {
    die "Invalid or missing input directory $inDir.";
} elsif (! -f "$inDir/6.1.1.20.fasta" || ! -f "$inDir/complete.genomes") {
    die "$inDir does not appear to contain a representative-genomes database.";
} else {
    print STDERR "Input database in $inDir.\n" if $debug;
}
# Determine the input style.
my ($p3, @filter, @cols, $roleCheck, $ih, $keyCol, $shrub);
my $fastaMode = $opt->fasta;
my $batchSize = $opt->batchsize;
# Get access to the databases.
if (! $opt->fasta) {
    print STDERR "Connecting to PATRIC.\n" if $debug;
    $p3 = P3DataAPI->new();
    print STDERR "Connecting to Shrub.\n" if $debug;
    $shrub = Shrub->new_for_script($opt);
    # Create the PATRIC filter and column clauses for genome queries.
    @filter = (['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']);
    @cols = qw(genome_id genome_name product aa_sequence);
    # Save the checksum for the seed role.
    $roleCheck = "WCzieTC/aZ6262l19bwqgw";
    # Open the input file.
    $ih = P3Utils::ih($opt);
    # Read the incoming headers to find the key column.
    (undef, $keyCol) = P3Utils::process_headers($ih, $opt);
} else {
    print STDERR "Opening FASTA input.\n" if $debug;
    my $fh = P3Utils::ih($opt);
    $ih = FastA->new($fh);
}
# Create the database from the input directory.
print STDERR "Creating database from $inDir.\n" if $debug;
my $repDB = RepGenomeDb->new_from_dir($inDir, unconnected => 1);
# Save the parameters.
my $K = $repDB->K();
my $minScore = $repDB->score();
print STDERR "Kmer size is $K and minimum similarity is $minScore.\n" if $debug;
# Write the output header.
P3Utils::print_cols([qw(genome_id genome_name rep_id rep_name score outlier)]);
# This will be a hash of bad genome IDs.
my %bad;
# Loop through the input.
my $done;
while (! $done) {
    # This will be a map of genomes to sequences we need to check, in the form [genome, name, null, protein].
    my %results;
    if ($fastaMode) {
        my $count = 0;
        while ($count < $batchSize && $ih->next) {
            my $genome = $ih->id;
            # If we are a feature ID, extract the genome ID from it.
            if ($genome =~ /^fig\|(\d+\.\d+)/) {
                $genome = $1;
            }
            # Store the sequence.
            $results{$genome} = [$genome, $ih->comment, '', $ih->left];
            $count++;
        }
        print STDERR "$count sequences read from FASTA.\n" if $debug;
    } else {
        # Here we are getting PATRIC/Shrub genome IDs.  This process is much more complex, since we have to sort out
        # the proteins we want from the extra junk returned by a typical query.
        my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
        my @genomes = map { $_->[0] } @$couplets;
        if (@genomes) {
            # First ask the Shrub for the names and identifying proteins of the un-represented genomes.
            # We keep the ones we don't find in here.
            my @rem;
            for my $genome (@genomes) {
                my @rows = $shrub->GetAll("Genome Genome2Feature Feature Feature2Function Function2Role AND Feature Feature2Protein Protein",
                              'Genome(id) = ? AND Feature2Function(security) = ? AND Function2Role(to-link) = ?',
                              [$genome,'0','PhenTrnaSyntAlph'], [qw(id name
                              Feature2Function(to-link) Protein(sequence))]);
                if (! @rows) {
                    push @rem, $genome;
                } else {
                    $stats->Add(genomeShrub => 1);
                    if (scalar(@rows) > 1) {
                        $stats->Add(badGenome => 1);
                        $bad{$genome} = 1;
                        $stats->Add(redundantProt => (scalar(@rows) - 1));
                    }
                    my $bestLen = 0;
                    for my $row (@rows) {
                        my (undef, $name, $function, $prot) = @$row;
                        my $protLen = length($prot);
                        if ($protLen > $bestLen) {
                            $results{$genome} = $row;
                            $bestLen = $protLen;
                        }
                    }
                }
            }
            # Ask PATRIC for the names and identifying proteins of the remaining genomes.
            my $resultList = P3Utils::get_data_keyed($p3, 'feature', \@filter, \@cols, \@rem, 'genome_id');
            # The resultList entries are in the form [$genome, $name, $function, $prot]. Get the longest
            # protein for each genome. We also track obviously bad genomes.
            my (%results);
            for my $result (@$resultList) {
                my ($genome, $name, $function, $prot) = @$result;
                # Check the protein.
                my $check = RoleParse::Checksum($function // '');
                if (! $prot) {
                    print STDERR "WARNING: $genome $name has no identifying protein.\n" if $debug;
                    $stats->Add(genomeNoProt => 1);
                    $stats->Add(badGenome => 1);
                    $bad{$genome} = 1;
                } elsif ($check ne $roleCheck) {
                    # Here the function matched but it is not really the same.
                    $stats->Add(funnyProt => 1);
                } else {
                    # Add the protein length to the result array.
                    my $protLen = length $prot;
                    push @$result, $protLen;
                    if (! exists $results{$genome}) {
                        $results{$genome} = $result;
                    } else {
                        $stats->Add(redundantProt => 1);
                        if (! $bad{$genome}) {
                            print STDERR "WARNING: $genome $name has a redundant identifying protein.\n" if $debug;
                            $bad{$genome} = 1;
                            $stats->Add(badGenome => 1);
                        }
                        if ($protLen > $results{$genome}[3]) {
                            # It's a better protein, so keep it.
                            $results{$genome} = $result;
                        }
                    }
                }
            }
        }
    }
    # Now loop through the genomes, checking for representatives.
    for my $genome (keys %results) {
        print STDERR "Checking $genome.\n" if $debug;
        my (undef, $name, undef, $prot) = @{$results{$genome}};
        my ($repID, $score) = $repDB->find_rep($prot);
        my $repName = '';
        $repID //= '';
        if ($repID) {
            $repName = $repDB->rep_object($repID)->name;
        }
        my $repFlag = '';
        if ($score < $minScore) {
            $stats->Add(genomeOutlier => 1);
            $repFlag = 'Y';
        }
        P3Utils::print_cols([$genome, $name, $repID, $repName, $score, $repFlag]);
    }
    # Check for end-of-file.
    if ($fastaMode) {
        $done = $ih->at_end();
    } else {
        $done = eof $ih;
    }
}
print STDERR "All done.\n" . $stats->Show() if $debug;


