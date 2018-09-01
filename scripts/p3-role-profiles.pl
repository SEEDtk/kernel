=head1 Role Profile Comparison

    p3-role-profiles.pl [options] genomeID <genomeList

This script compares the key roles in a specified input genome with the roles in a list of theoretically close genomes
specified in the standard input file.  For each close genome, the output file will list the roles it has that the
specified genome does not have and the roles the specified genome has that the close genome does not have.

=head2 Parameters

The positional parameter is the ID of the genome of interest.

The standard input can be overridden using the options in L<P3Utils/ih_options>. The column specifying the genome ID
can be selected using the options in L<P3Utils/col_options>.  The delimiter to use between role IDs in the output can
be specified with L<P3Utils/delim_options>.  The role specification files can be configured using the options if
L<EvalCon/role_options>.

The following additional command-line options are supported.

=over 4

=item verbose

Write status information to the standard error output.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use EvalCon;
use GEO;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('genomeID', P3Utils::col_options(), P3Utils::ih_options(),
        P3Utils::delim_options(), EvalCon::role_options(),
        ['verbose|debug|v', 'show status messages on STDERR'],
        );
my $stats = Stats->new();
# Get the list delimiter.
my $delim = P3Utils::delim($opt);
# Get the debug flag.
my $debug = $opt->verbose;
my $logH = ($debug ? \*STDERR : undef);
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Read the role definitions.
my $evalCon = EvalCon->new_for_script($opt, $logH);
my ($nMap, $cMap) = $evalCon->roleHashes;
my @roleList = sort @{$evalCon->rolesToUse};
my %rolesToUseH = map { $_ => 1 } @roleList;
# Compute the options for GEO creation.
my %geoOptions = (roleHashes => [$nMap, $cMap], stats => $stats, detail => 0, logH => $logH, rolesToUse => \%rolesToUseH);
# Get the input genome ID and extract its roles.
my $roleHash;
my ($genomeID) = @ARGV;
if (! $genomeID) {
    die "No target genome ID specified."
} else {
    my $gHash = GEO->CreateFromPatric([$genomeID], %geoOptions);
    my $geo = $gHash->{$genomeID};
    if (! $geo) {
        die "Genome $genomeID not found in PATRIC.";
    } else {
        $roleHash = $geo->roleCounts;
    }
}
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders, $keyCol) = P3Utils::process_headers($ih, $opt);
# Form the full header set and write it out.
if (! $opt->nohead) {
    P3Utils::print_cols([@$outHeaders, 'new_roles', 'missing_roles']);
}
# This hash will be used for the role report. The role report has one row for each role
# and one column for each genome.
my %roleGenomes;
# Loop through the input.
while (! eof $ih) {
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # Get the GEOs for these genomes.
    my @gList = map { $_->[0] } @$couplets;
    my $gHash = GEO->CreateFromPatric(\@gList, %geoOptions);
    # Loop through the couplets, doing the comparison.
    for my $couplet (@$couplets) {
        my ($genome, $line) = @$couplet;
        my $geo = $gHash->{$genome};
        if ($geo) {
            # Here we found the genome and we can compare the roles.
            my (@missing, @new);
            my $roleH = $geo->roleCounts;
            for my $role (@roleList) {
                my $new1 = $roleH->{$role};
                my $old1 = $roleHash->{$role};
                if (! $new1 && $old1) {
                    push @missing, $role;
                    $stats->Add(missingRole => 1);
                    $roleGenomes{$role}{$genome} = '-';
                } elsif ($new1 && ! $old1) {
                    push @new, $role;
                    $stats->Add(newRole => 1);
                    $roleGenomes{$role}{$genome} = '+';
                }
            }
            P3Utils::print_cols([@$line, scalar(@new), scalar(@missing)]);
            for my $type ([new => @new], [missing => @missing]) {
                my ($label, @roles) = @$type;
                for my $role (@roles) {
                    P3Utils::print_cols(['', $label, $role, $nMap->{$role}]);
                }
            }
        }
    }
    # Now create the role matrix.
    P3Utils::print_cols([]);
    P3Utils::print_cols(['role', @gList]);
    for my $role (sort keys %roleGenomes) {
        my $gData = $roleGenomes{$role};
        P3Utils::print_cols([$nMap->{$role}, map { $gData->{$_} // '' } @gList]);
    }
    P3Utils::print_cols([]);
    P3Utils::print_cols([]);
}
if ($debug) {
    print STDERR "All done.\n" . $stats->Show();
}
