=head1 Find Best Representative Genomes

    rep_update.pl [options] sourceDir targetDir

This rather esoteric script takes the representative-genomes from a source directory and finds which genomes best represent them
in a target directory. The genomes are added to the target directory's represented-genomes database.

=head2 Parameters

The positional parameters are the names of the source directory and the target directory. Both are representative-genome
directories as described in L<RepGenomeDb/new_from_dir>. The source directory does not need a C<rep_db.tbl> file. If
the C<--clear> option is specified, the C<rep_db.tbl> in the target directory will be deleted before loading.

The command-line options are as follows.

=over 4

=item clear

Erase the represented-genomes file in the target directory before starting. Only the new genomes will be in the output.

=back

=cut

use strict;
use P3Utils;
use RepGenomeDb;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('sourceDir targetDir',
        ['clear', 'erase old represented-genomes information before processing']
        );
my $stats = Stats->new();
# Get the input directories.
my ($sourceDir, $targetDir) = @ARGV;
if (! $sourceDir) {
    die "No source directory specified.";
} elsif (! -d $sourceDir) {
    die "Source directory $sourceDir not found or invalid.";
} elsif (! $targetDir) {
    die "No target directory specified.";
} elsif (! -d $targetDir) {
    die "Target directory $targetDir not found or invalid.";
}
# Load the source directory.
print "Loading source directory $sourceDir.\n";
my $sourceDb = RepGenomeDb->new_from_dir($sourceDir, unconnected => 1, verbose => 1);
# Load the target directory.
my $clear = ($opt->clear ? 1 : 0);
print "Loading target directory $targetDir.\n";
my $targetDb = RepGenomeDb->new_from_dir($targetDir, unconnected => $clear, verbose => 1);
# Get the list of source genomes.
my $sGenomes = $sourceDb->rep_list();
# Get a list of the genomes from the target list.
my %tGenomes = map { $_ => 1 } @{$targetDb->rep_list()};
# Get the match score for the target database.
my $score = $targetDb->score();
print "Target database score is $score.\n";
# Loop through the source genomes.
my $count = 0;
my $found = 0;
for my $genome (@$sGenomes) {
    # Insure this genome is not already a target.
    my $skip;
    if ($tGenomes{$genome}) {
        $stats->Add(genomeRepresentative => 1);
        $skip = 1;
    } elsif (! $clear) {
        # Insure it is not already represented.
        my ($fID, $fScore) = $targetDb->check_rep($genome);
        if ($fID) {
            $stats->Add(genomeAlreadyRepresented => 1);
            $skip = 1;
        }
    }
    if (! $skip) {
        # Get the genome's protein.
        my $gObject = $sourceDb->rep_object($genome);
        my $prot = $gObject->prot();
        # Find the best match.
        my ($fID, $fScore) = $targetDb->find_rep($prot, $score);
        if ($fScore < $score) {
            my $name = $gObject->name;
            print "$genome ($name) has no representatives.\n";
            $stats->Add(genomeNotFound => 1);
        } else {
            $targetDb->Connect($fID, $genome, $fScore);
            $stats->Add(genomeFound => 1);
            $found++;
        }
    }
    $count++;
    print "$count genomes processed, $found representatives found.\n" if $count % 100 == 0;
}
print "Saving $targetDir.\n";
$targetDb->Save($targetDir);
print "All done.\n" . $stats->Show();
