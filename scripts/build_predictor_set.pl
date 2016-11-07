#!/usr/bin/env perl

use IPC::Run3;

use strict;
use warnings;
use Data::Dumper;
use Carp;

use IO::File;
use File::Spec;
use File::Path qw(make_path remove_tree);
require Proc::ParallelLoop;

use SeedUtils;

my ($data_dir, $output_dir, $fraction, $classifier) = @ARGV;

if (!-d $output_dir) {
    &run_safe(['cp', '-pR', $data_dir, $output_dir], \undef, \undef, \*STDERR);
}
else {
    die "Output-directory '$output_dir' already exists";
}

my @funcs = map { m/^\d+\t(\S+)/ ? $1 : () } SeedUtils::file_read("$data_dir/col.h");
chomp @funcs;

my %func_to_col = map { chomp; m/^(\d+)\t(\S+)/ ? ($2 => $1) : () } &SeedUtils::file_read("$data_dir/col.h");

Proc::ParallelLoop::pareach(\@funcs, \&worker, {Max_Workers => 32});

sub worker {
    my ($func) = @_;

    my $col = $func_to_col{$func};
    unless (defined($col)) {
        print STDERR "Could not find column for func='$func' -- skipping\n";
        next;
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
