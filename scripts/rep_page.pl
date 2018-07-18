=head1 Create Representative-Genome Web Page

    rep_page.pl [options] repDir outFile

This script produces a web page for examining the representative genomes and their members.

=head2 Parameters

The positional parameters are the name of the representative-genome directory and the name of the output web page file.
The representative-genome directory is described in detail in L<RepGenomeDb/new_from_dir>. The command-line options are
as follows.

=over 4

=item template

The name of the web page template file. The default is C<reps.tt> in the L<BinningReports> template directory.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use RepGenomeDb;
use Template;
use Stats;
use RoleParse;
use URI::Escape;

# Get the command-line options.
my $opt = P3Utils::script_opts('repDir outFile',
        ['template=s', 'web page template file', { default => "$FIG_Config::mod_base/RASTtk/lib/BinningReports/reps.tt" }],
        );
my $stats = Stats->new();
# Get the input directory and the output file.
my ($repDir, $outFile) = @ARGV;
if (! $repDir) {
    die "No input directory specified.";
} elsif (! -d $repDir) {
    die "Input directory $repDir not found or invalid.";
} elsif (! $outFile) {
    die "No output file specified.";
}
# Verify that we can open the output file.
open(my $oh, ">$outFile") || die "Could not open output file: $!";
# Get access to PATRIC.
print "Connecting to PATRIC.\n";
my $p3 = P3DataAPI->new();
# Load the rep-genomes database.
my $repDB = RepGenomeDb->new_from_dir($repDir, verbose => 1);
# The template variables will be built in here.
my %vars = ( kSize => $repDB->K, score => $repDB->score );
# Get the full list of representative genomes and store the count.
my $genomesL = $repDB->rep_list;
$vars{reps} = scalar @$genomesL;
# This hash will contain a descriptor for each genome.
my %genomes;
# Get key information about each of the genomes from PATRIC.
print "Retrieving genome statistics.\n";
my $genomeDataL = P3Utils::get_data_keyed($p3, genome => [], ['genome_id', 'genome_name', 'contigs', 'genome_length', 'kingdom'], $genomesL);
for my $genomeDatum (@$genomeDataL) {
    my ($id, $name, $contigs, $dna_bp, $kingdom) = @$genomeDatum;
    $genomes{$id} = { id => $id, name => $name, contigs => $contigs, dna_bp => $dna_bp, kingdom => $kingdom, seed => '&nbsp;' };
    $stats->Add(genomeDatum => 1);
}
# Get the seed protein IDs from PATRIC. We know there will be only one, but we have to check to insure we don't get similar function
# names mixed in.
print "Retrieving seed protein IDs.\n";
$genomeDataL = P3Utils::get_data_keyed($p3, feature => [['eq', 'product', 'Phenylalanyl-tRNA synthetase alpha chain']],
        ['genome_id', 'patric_id', 'product'], $genomesL, 'genome_id');
for my $genomeDatum (@$genomeDataL) {
    my ($id, $fid, $function) = @$genomeDatum;
    $stats->Add(featureDatum => 1);
    my $checksum = RoleParse::Checksum($function);
    if ($checksum eq 'WCzieTC/aZ6262l19bwqgw') {
        # Here we need to create the feature link.
        my $url = "https://www.patricbrc.org/view/Feature/" . uri_escape($fid);
        $genomes{$id}{seed} = qq(<a href="$url" target="_blank">$fid</a>);
        $stats->Add(seedProteinFound => 1);
    } else {
        $stats->Add(funnySeed => 1);
    }
}
# Now loop through the genomes, creating the member lists. We will output the fully-formed genome descriptors in this list.
# Our main task here is to create the membership string.
print "Assembling genome descriptors.\n";
my $mTotal = 0;
for my $genome (@$genomesL) {
    $stats->Add(genomeProcessed => 1);
    my $repObject = $repDB->rep_object($genome);
    my $descriptor = $genomes{$genome};
    # Get the member list. What comes back is [ID, score] 2-tuples. We sort them and extract the IDs.
    my $repList = $repObject->rep_list();
    my @members = map { $_->[0] } sort { $b->[1] <=> $a->[1] } @$repList;
    my $mCount = scalar @members;
    if (! $mCount) {
        $descriptor->{membership} = '&nbsp;';
    } elsif ($mCount == 1) {
        # Here we have 1 member.
        $descriptor->{membership} = qq(<a href="https://patricbrc.org/view/Genome/$genome" target="_blank">1 member</a>);
    } elsif ($mCount < 700) {
        # Here we have a bunch of members.
        my $genomeIDs = join(",", @members);
        $descriptor->{membership} = qq(<a href="http://patricbrc.org/view/GenomeList?in(genome_id,($genomeIDs))" target="_blank">$mCount members</a>);
    } else {
        # Here the list is too big, so we have to break it into pieces.
        my @links;
        my $pageI = 1;
        for (my $i = 0; $i < $mCount; $i += 700) {
            my $j = $i + 699;
            $j = $mCount - 1 if $j >= $mCount;
            my $genomeIDs = join(",", @members[$i .. $j]);
            push @links, qq(<a href="http://patricbrc.org/view/GenomeList?in(genome_id,($genomeIDs))" target="_blank">$pageI</a>);
            $pageI++;
        }
        $descriptor->{membership} = join(' ', "$mCount members:", @links);
    }
    # Save the member count.
    $descriptor->{mCount} = $mCount;
    $mTotal += $mCount;
}
# Store the member total and the genome list in the variable hash.
$vars{genomes} = [ map { $genomes{$_} } sort { $genomes{$b}{mCount} <=> $genomes{$a}{mCount} } keys %genomes ];
$vars{mTotal} = $mTotal;
# We have now compiled the information we need for the report. Create the template engine.
print "Constructing web page in $outFile.\n";
my $templateEngine = Template->new(ABSOLUTE => 1);
# Allocate the result variable.
my $html;
# Create the web page.
my $template = $opt->template;
print "Tempate file is $template\n.";
$templateEngine->process($opt->template, \%vars, \$html) || die $templateEngine->error();
# Output the web page.
print $oh $html;
# All done.
print "All done.\n" . $stats->Show();