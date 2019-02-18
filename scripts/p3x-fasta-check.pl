=head1 Verify Genomes in a FASTA File

    p3x-fasta-check.pl [options]

This is a simple program that runs through a FASTA file, parsing the genome ID out of the ID in each record, and queries the
PATRIC database to make sure the genome is valid.

=head2 Parameters

There are no positional parameters.

The standard input contains the FASTA file and can be overridden using the options in L<P3Utils/ih_options>.

The command-line options are as follows.

=over 4

=item verbose

Write status messages to STDERR.

=item corrected

If specified, a file name.  A corrected FASTA file will be written to this file name.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use Stats;
use FastA;

# Get the command-line options.
my $opt = P3Utils::script_opts('', P3Utils::ih_options(),
            ['verbose|debug|v', 'write status messages to STDERR'],
            ['corrected=s', 'if specified, the name of a file to which the corrected FASTA data should be written']
        );
my $stats = Stats->new();
# Check for a correction file.
my $oh;
if ($opt->corrected) {
    open($oh, '>', $opt->corrected) || die "Could not open corrections file: $!";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Open the input file.
my $ih = P3Utils::ih($opt);
# Check the options.
my $debug = $opt->verbose;
# Write the output headers.
P3Utils::print_cols(['genome', 'comment', 'found', 'name']);
# Loop through the FASTA file.
my %genomes;
my $count = 0;
print STDERR "Looping through FASTA file.\n" if $debug;
my $faH = FastA->new($ih);
while ($faH->next()) {
    $genomes{$faH->id} = [$faH->comment, $faH->left];
    $stats->Add(genomeIn => 1);
    $count++;
    if ($count >= 100) {
        ProcessBatch(\%genomes);
        %genomes = ();
        $count = 0;
    }
}
if ($count > 0) {
    ProcessBatch(\%genomes);
}
print STDERR "All done.\n" . $stats->Show() if $debug;

# Process a batch of genome IDs.
sub ProcessBatch {
    my ($genomeH) = @_;
    my @ids = sort keys %$genomeH;
    print STDERR "Fetching genome names from PATRIC.\n" if $debug;
    my $genomeData = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name'], \@ids);
    print STDERR scalar(@$genomeData) . " names found of " . scalar(@ids) . "\n" if $debug;
    # Form the genomes found into a hash.
    my %found;
    for my $genomeDatum (@$genomeData) {
        my ($genome, $name) = @$genomeDatum;
        $found{$genome} = $name;
        $stats->Add(genomeFound => 1);
    }
    # Verify we found them all.
    for my $id (@ids) {
        my ($found, $name, $comment, $seq) = ('', '', @{$genomeH->{$id}});
        if ($found{$id}) {
            $found = 1;
            $name = $found{$id};
            if ($oh) {
                print $oh ">$id $comment\n$seq\n";
            }
        } else {
            $stats->Add(genomeNotFound => 1);
        }
        P3Utils::print_cols([$id, $comment, $found, $name]);
    }
    $stats->Add(batchProcessed => 1);
}