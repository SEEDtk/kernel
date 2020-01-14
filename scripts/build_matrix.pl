#!/usr/bin/env perl

=head1 Create Probdir: Function Predictors Step 2

    build_matrix raw.table probDir

This script takes the output from L<build_role_tables.pl> or L<p3-build-role-tables.pl> and creates a function matrix. Each row of the
matrix represents a genome. Each column of the matrix represents a role. The values in the matrix indicate
how many times the role occurs in the genome. The resulting matrix is encoded as a SciKit probdir. In other
words, there is a C<row.h> of genome IDs, an C<col.h> of role IDs, and an C<X> file containing the actual matrix.

This is the second step of the process for building function predictors. The third is B<build_LDA_predictors.py>.
The final step is L<build_roles_to_use.pl>, after which this script is used again if the predictors have not stabilized.

=head2 Parameters

The positional parameters are the name of the C<raw.table> file and the name of the output directory.

The output directory cannot exist unless the C<--clear> option is specified. If it does exist, the C<roles.to.use>
and C<roles.in.subsystems> files will be taken from it if they exist. Otherwise they will default to the global files.

=cut

use strict;
use warnings;
use Data::Dumper;
use ScriptUtils;
use SeedUtils;
use File::Copy::Recursive;
use EvalCon;

$| = 1;
my %funcs;
my %genomes;
my %counts;

my $opt = ScriptUtils::Opts('raw.table probDir',
        ['clear', 'overwrite previous results']);

my ($table_file, $probDir) = @ARGV;
if (! $table_file) {
    die "No raw.table specified.";
} elsif (! -s $table_file) {
    die "$table_file is missing or empty.";
} elsif (! $probDir) {
    die "No output directory specified.";
} elsif (-f $probDir) {
    die "$probDir is not a directory.";
}


if (!-d $probDir) {
    File::Copy::Recursive::pathmk($probDir) or die "Could not create probDir '$probDir': $!";
} elsif ($opt->clear) {
    print "Clearing $probDir.\n";
    if (-d "$probDir/Predictors") {
        File::Copy::Recursive::pathrmdir("$probDir/Predictors") || die "Could not clear Predictors directory: $!";
    }
    my @files = qw(col.h roles.mapped roles.not_mapped row.h X);
    for my $file (@files) {
        if (-f "$probDir/$file") {
            unlink "$probDir/$file";
        }
    }
} else {
    die "ERROR: probDir '$probDir' already exists";
}

# Create the eval object.
my $eval = EvalCon->new(predictors => $probDir, logH => \*STDOUT);
# Process the input and produce the output.
print "Processing $table_file.\n";
$eval->BuildFileMatrix($table_file, $probDir);
print "All done.\n" . $eval->stats->Show();

