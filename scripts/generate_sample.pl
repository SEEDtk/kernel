use strict;
use Data::Dumper;
use List::Util qw(shuffle);
use Shrub::Contigs;
use Shrub;
use ScriptUtils;
use File::Slurp;

=head1 Generate Community Pipeline Sample

    generate_sample [ options ] coverage genomes >fasta

Generate a FASTA file of contigs as test data for L<community_pipeline.pl>. A random sample of genomes will be selected from
the Shrub database. Alternatively, an input file of genome IDs can be specified so that specific genomes will be selected.
The FASTA output will contain contigs created from the sample genomes. In addition, a blacklist file will be produced
listing the genomes used, and a bin file will be produced that relates each contig to the genome from which it came. These
last two files are used as input to community_pipeline:  the blacklist file to C<--blacklist> and the bin file to
C<--expected>.

=head2 Parameters

There are two positional parameters-- a maximum coverage and a genome specification. The maximum coverage is
a number that must be greater than 10. The genome specification can be a number-- in which case that number of genomes will
be randomly selected-- or the name of a file containing genome IDs to use. The file must be tab-delimited, each record
containing (0) a genome ID and (1) a genome name.

The command-line options are those found in L<Shrub/script_options> plus the following.

=over 4

=item blacklist

If specified, the name of the file to contain the list of IDs for the genomes used. This will be a tab-delimited file,
each record containing (0) the genome ID, (1) the genome name, and (2) the assigned coverage for that genome's contigs.

=item expected

If specified, the name of the file to contain the bin information. This will be a two-column tab-delimited file, each
record containing (0) a contig ID and (1) the ID of the genome from which it came.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('abundance genomes', Shrub::script_options(),
        ['blacklist|b=s', 'output file for list of genomes used'],
        ['expected|e=s', 'output file for expected result bins']
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);

my ($Na, $genomes) = @ARGV;
if (! $Na || $Na < 10) {
    die "MaxAbundance must be 10 or more";
}
my @genomes;
if (! $genomes) {
    die "No genome parameter specified.";
} elsif ($genomes =~ /^\d+$/) {
    @genomes = shuffle $shrub->GetAll('Genome', '', [], 'id name');
    $#genomes = $genomes-1;
} else {
    @genomes = map { [ split /\t/, $_ ] } File::Slurp::read_file($genomes, { chomp => 1 });
}
my ($bh, $eh);
if ($opt->blacklist) {
    open($bh, ">", $opt->blacklist) || die "Could not open the blacklist output file: $!";
}
if ($opt->expected) {
    open($eh, ">", $opt->expected) || die "Could not open the bin output file (expected): $!";
}
my $id = 1;
foreach my $gData (@genomes)
{
    my ($g, $gname) = @$gData;
    my $contigO = Shrub::Contigs->new($shrub,$g);
    my @triples = $contigO->tuples;
    my $real_ln = 0;
    foreach $_ (@triples)
    {
        $real_ln += length($_->[2]);
    }
    my $cov = 10 + int(rand() * ($Na - 10));
    Print($bh, join("\t", $g, $gname, $cov), "\n");
    my $out_ln = 0;
    my $copies = 0.5 + rand();
    while ($out_ln < ($copies * $real_ln))
    {
        my $desired = int(rand() * 10000) + 200;
        my $start = int(rand() * scalar @triples);
        my $got = 0;
        while (! $got)
        {
            my $len_of_contig = length($triples[$start]->[2]);
            if ($desired > $len_of_contig)
            {
                $start++;
                if ($start == @triples) { $start = 0 }
            }
            else
            {
                my $first = int(rand() * ($len_of_contig - $desired));
                my $seq = substr($triples[$start]->[2],$first,$desired);
                my $contigID = join("_", $g, $id++, "length", $desired, "cov", $cov);
                print ">", $contigID, "\n";
                print $seq,"\n";
                Print($eh, join("\t", $contigID, $g), "\n");
                $out_ln += $desired;
                $got = 1;
            }
        }
    }
}

# Print to a file handle if it is defined.
sub Print {
    my ($h,@text) = @_;
    if (defined $h) {
        print $h @text;
    }
}