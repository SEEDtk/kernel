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

=item minRoles

Minimum number of acceptable singly-occurring roles. The default is C<1000>.

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

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outDir', P3Utils::col_options(), P3Utils::ih_options(),
        ['roleFile|rolefile|R=s', 'role definition file', { default => "$FIG_Config::global/roles.in.subsystems"} ],
        ['minRoles|minroles|m=i', 'minimum number of acceptable roles per genome', { default => 1000 }]
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
open(my $oh, ">$outDir/genomes.tbl") || die "Could not open genomes.tbl: $!";
P3Utils::print_cols(['genome_id', 'lineage', 'roles'], oh => $oh);
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    my $gHash = GEO->CreateFromPatric([map { $_->[0] } @$couplets], %geoOptions);
    print scalar(keys %$gHash) . " genomes found in batch.\n";
    my $kept = 0;
    # Loop through the individual genomes.
    for my $genome (keys %$gHash) {
        my $geo = $gHash->{$genome};
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
    }
    print "$kept genomes output.\n";
}
print "All done.\n" . $stats->Show();