#!/usr/bin/env perl

=head1 Create Predictors: Function Predictors Step 3

    build_predictor_set [options] probDir

This is the third step in creating function predictors, after L<build_matrix.pl>. It takes as input the probdir
from that script and creates predictors for the individual roles in the subdirectory C<Predictors>. After this,
use L<build_roles_to_use.pl>.

=head2 Parameters

The position parameter is the name of the function predictor directory, containing the initial files created by
L<build_matrix.pl>.

The command-line options are

=over 4

=item clear

If specified, the C<Predictors> directory will be erased if it currently exists; otherwise, an existing output directory
will cause an error.

=item fraction

The fraction of matrix entries to use for testing instead of training. The default is C<0.2>, indicating one in
five.

=item classifier

The type of classifier to use. The default is C<RandomForestClassifier>.

=back

=cut

use IPC::Run3;

use strict;
use warnings;
use Data::Dumper;
use Carp;

use IO::File;
use File::Spec;
use File::Path qw(make_path remove_tree);
require Proc::ParallelLoop;
use ScriptUtils;
use SeedUtils;
use File::Copy::Recursive;

$| = 1;
my $opt = ScriptUtils::Opts('probDir',
        ['fraction|f=f', 'fraction of entries to keep in random samples', { default => 0.2 }],
        ['classifier|c=s', 'type of classifier to use', { default => 'RandomForestClassifier'}],
        ['clear', 'erase output directory before starting']);
my ($data_dir) = @ARGV;
my $fraction = $opt->fraction;
my $classifier = $opt->classifier;

if (! -d $data_dir) {
    die "Invalid data directory.";
} elsif (! -s "$data_dir/X") {
    die "build_matrix has apparently not been run on $data_dir.";
}
my $output_dir = "$data_dir/Predictors";
if (!-d $output_dir) {
    mkdir $output_dir;
} elsif ($opt->clear) {
    print STDERR "Clearing output directory.\n";
    File::Copy::Recursive::pathempty($output_dir);
}
# Get the current directories so we can skip them.
print "Checking output directory.\n";
opendir(my $dh, $output_dir) || die "Could not open $output_dir: $!";
my %already = map { $_ => 1 } grep { -s "$output_dir/$_/Classifiers/$classifier/accuracy" } readdir $dh;
closedir $dh;
print STDERR scalar(keys %already) . " functions already processed.\n";
my @funcs = grep { ! $already{$_} } map { m/^\d+\t(\S+)/ ? $1 : () } SeedUtils::file_read("$data_dir/col.h");
chomp @funcs;
print STDERR scalar(@funcs) . " functions in table.\n";

my %func_to_col = map { chomp; m/^(\d+)\t(\S+)/ ? ($2 => $1) : () } &SeedUtils::file_read("$data_dir/col.h");

Proc::ParallelLoop::pareach(\@funcs, \&worker, {Max_Workers => 32});

sub worker {
    my ($func) = @_;
    my $col = $func_to_col{$func};
    unless (defined($col)) {
        print STDERR "Could not find column for func='$func' -- skipping\n";
        next;
    } else {
        print STDERR "Processing $col: $func.\n";
    }

    my $probDir = "$output_dir/$func";
    &run_safe(["partition_on_column", $data_dir, $probDir, $col],
              \undef, \undef, \*STDERR
        );

    make_path("$probDir/Classifiers/$classifier");
    my $accuracy = "$probDir/Classifiers/$classifier/accuracy";
    &run_safe(["classifier_cross_validation",
               "--probDir=$probDir",
               "--fraction=$fraction",
               "--classifier=$classifier",
               "--output=$accuracy",
               "--error=$probDir/Classifiers/$classifier/cross-validation.err"
              ],
              \undef, \undef, \*STDERR
    );

    my $train_err_fh = new IO::File ">> $probDir/Classifiers/$classifier/final_training.log";
    &run_safe(["train_classifier", $probDir, $classifier],
        \undef, \undef, $train_err_fh
        );
}
print STDERR "All done.\n";
exit(0);



sub run_safe {
    my ( $args, $in_fh, $out_fh, $err_fh ) = @_;
#   print $err_fh Dumper($in_fh, $out_fh, $err_fh);

    if (my $rc = run3( $args, $in_fh, $out_fh, $err_fh )) {
        return $rc;
    }
    else {
        if ($? == -1) {
            print $err_fh "failed to execute: $!\n";
            confess("aborting");
        }
        elsif ($? & 127) {
            print $err_fh ("child died with signal %d, %s coredump\n",
                           ($? & 127),  ($? & 128) ? 'with' : 'without'
                );
            confess("aborting");
        }
    }
}
