=head1 Compute Universal Roles and Lineages

    p3-genome-unis.pl [options] outDir

This script takes a file of genome IDs as input and outputs the universal roles and taxonomic lineage of each. The standard input must contain
genome IDs in the key column. The output directory will be tab-delimited. The first column will contain a genome ID, the second a comma-delimited
list consisting of the taxonomic lineage (from largest group to smallest), and the third a comma-delimited list of singly-occurring role IDs.

If a genome has too few singly-occurring roles it will be skipped.

=head2 Parameters

The positional parameter is the name of the output directory. The output file will be named C<genomes.tbl>.

The standard input can be overridden using the options in L<P3Utils/ih_options>. The key column and batch size can be specified using the
options in L<P3Utils/col_options>.

Additional command-line options are the following.

=over 4

=item resume

Resume output after the specified genome. Use this to restart after a failure.

=item minRoles

Minimum number of acceptable singly-occurring roles. The default is C<800>.

=item roleFile

The C<roles.in.subsystems> file containing the role IDs, checksums, and names for the stable roles. The default is the
one in the SEEDtk global directory.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use GEO;
use EvalCon;
use Stats;
use File::Copy::Recursive;
use Time::HiRes;
use Math::Round;
use IO::Handle;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir', P3Utils::col_options(), P3Utils::ih_options(),
        ['roleFile|rolefile|R=s', 'role definition file', { default => "$FIG_Config::global/roles.in.subsystems"} ],
        ['minRoles|minroles|m=i', 'minimum number of acceptable roles per genome', { default => 800 }],
        ['resume=s', 'resume after the specified genome']
        );
my $stats = Stats->new();
# Get the output directory.
my ($outDir) = @ARGV;
if (! $outDir) {
    die "Could not open output directory.";
} elsif (! -d $outDir) {
    File::Copy::Recursive::pathmk($outDir) || die "Could not create output directory $outDir: $!";
}
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Read the role file.
my $roleFile = $opt->rolefile;
print "Reading role file $roleFile.\n";
my ($nMap, $cMap) = EvalCon::LoadRoleHashes($roleFile, $stats);
# Create the GEO options.
my %geoOptions = (detail => 0, roleHashes => [$nMap, $cMap], p3 => $p3, stats => $stats, logH => \*STDOUT);
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Open the output file.
my $mode = ($opt->resume ? '>>' : '>');
open(my $oh, $mode, "$outDir/genomes.tbl") || die "Could not open genomes.tbl: $!";
$oh->autoflush(1);
my $lastGenome;
if (! $opt->resume) {
    P3Utils::print_cols(['genome_id', 'lineage', 'roles'], oh => $oh);
} else {
    # Here we have to resume in the middle.
    my $found;
    my $count = 0;
    $lastGenome = $opt->resume;
    print "Searching for $lastGenome.\n";
    while (! eof $ih && ! $found) {
        my $line = <$ih>;
        my @fields = P3Utils::get_fields($line);
        $count++;
        if ($fields[$keyCol] eq $lastGenome) {
            $found = 1;
        }
    }
    print "$count lines skipped.\n";
}
# Loop through the input.
eval {
    my $start0 = time;
    my ($count, $ocount) = (0,0);
    while (! eof $ih) {
        my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
        my $gHash = GEO->CreateFromPatric([map { $_->[0] } @$couplets], %geoOptions);
        print scalar(keys %$gHash) . " genomes found in batch.\n";
        my $kept = 0;
        # Loop through the individual genomes.
        for my $genome (keys %$gHash) {
            my $geo = $gHash->{$genome};
            $count++;
            # Get the lineage.
            my @lineage = @{$geo->lineage};
            # Compute the universal roles.
            my $roleH = $geo->roleCounts;
            my @unis = sort grep { $roleH->{$_} == 1 } keys %$roleH;
            $stats->Add(genomesProcessed => 1);
            my $uniCount = scalar @unis;
            if ($uniCount < $opt->minroles) {
                $stats->Add(genomeRejected => 1);
                print "$genome has only $uniCount roles-- skipped.\n";
            } else {
                $stats->Add(unisFound => scalar @unis);
                # Write the output.
                P3Utils::print_cols([$genome, join(",", @lineage), join(",", @unis)], oh => $oh);
                $kept++;
            }
            $lastGenome = $genome;
        }
        $ocount += $kept;
        print "$kept genomes output in this batch. $ocount output of $count processed at " . Math::Round::nearest(0.01, (time - $start0) / $count) . " seconds/genome.\n";
    }
};
if ($@) {
    if (! $lastGenome) {
        print "ERROR before any output.\n";
    } else {
        print "ERROR after $lastGenome.\n";
    }
    print "$@\n";
}
print "All done.\n" . $stats->Show();