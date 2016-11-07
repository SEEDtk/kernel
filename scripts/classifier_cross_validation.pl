#!/usr/bin/env perl

use IPC::Run3;

use strict;
use warnings;
use Data::Dumper;
use Carp;

use IO::File;
use File::Spec;
use File::Basename qw(basename);
use File::Path qw( make_path remove_tree );
use FIG_Config;

my $usage;

my $probDir;
my $fraction   = 0.2;
my $num_sweeps = 10;
my $classifier = 'RandomForestClassifier';
my $output;
my $error;
my $help;

use Getopt::Long;
my $rc = GetOptions(
    'probDir=s'    => \$probDir,
    'fraction=f'   => \$fraction,
    'sweeps=i'     => \$num_sweeps,
    'classifier=s' => \$classifier,
    'output=s'     => \$output,
    'error=s'      => \$error,
    );
print STDERR qq(\nrc=$rc\n\n) if $ENV{VERBOSE};

if (!$rc || $help) {
    warn qq(\n   usage: $usage\n\n);
    exit(0);
}

unless ($probDir && $fraction && $classifier) {
    if (@ARGV == 3) {
        ($probDir, $fraction, $classifier) = @ARGV;
    }
    else {
        die "Missing arguments";
    }
}


my $out_fh;
if ($output) {
    open($out_fh, '>', $output)
        or die "Could not write-open '$output'";
}
else {
    $out_fh = \*STDOUT;
}

my $err_fh;
if ($error)  {
    open($err_fh, '>', $error)
        or die "Could not write-open '$error'";
}
else {
    $err_fh  = \*STDERR;
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Define scratch directory names...
#-----------------------------------------------------------------------
my $trainDir = $FIG_Config::temp . '/tmp.' . basename($probDir) . q(.train.) . $$;
my $testDir  = $FIG_Config::temp . '/tmp.' . basename($probDir) . q(.test.) . $$;
my $result_file  = $FIG_Config::temp . q(/tmp.apply.result.) . $$;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Define training logfile...
#-----------------------------------------------------------------------
make_path("$probDir/Classifiers/$classifier");

my $train_err_fh = new IO::File;
my $training_logfile = "$probDir/Classifiers/$classifier/cross-validation_training.log";
$train_err_fh->open("> $training_logfile")
    or die "Could not write-open '$training_logfile'";

my @acc;
for (my $sweep=1; $sweep <= $num_sweeps; ++$sweep) {
    print $train_err_fh "Beginning sweep=$sweep\n";

    remove_tree($trainDir, $testDir);

    print $err_fh "Doing split\n" if $ENV{VERBOSE};
    &run_safe([ "split_probdir", $probDir, $trainDir, $testDir, $fraction ],
              \undef, \undef, $err_fh
        ) || die "Split failed on sweep=$sweep: $?, $!";

    print $err_fh "Doing train\n" if $ENV{VERBOSE};
    &run_safe([ "train_classifier", $trainDir, $classifier ],
              \undef, \undef, $train_err_fh
        ) || die "Training failed on sweep=$sweep: $?, $!";

    print $err_fh "Doing apply\n" if $ENV{VERBOSE};
    my $apply_out_fh = new IO::File;
    $apply_out_fh->open( qq(> $result_file) );
    &run_safe([ "apply_classifier", $testDir, $trainDir, $classifier ],
              \undef, $apply_out_fh, $train_err_fh
        ) || die "Application failed on sweep=$sweep: $?, $!";
    $apply_out_fh->close();
    $apply_out_fh->open( qq(< $result_file) );
    my @result = <$apply_out_fh>;
    $apply_out_fh->close();

    my $num_right = 0;
    my $num_instances = @result;
    map { chomp;
          my ($idx, $true, $pred) = split /\t/;
          if ($true == $pred) { ++$num_right; }
    } @result;
    my $acc = sprintf("%.2f", 100.0 * $num_right / $num_instances );

    push @acc, $acc;
    print $train_err_fh "$acc\n\n";
}
@acc = map { sprintf('%.2f', $_) } sort { $a <=> $b } @acc;


my $mean = 0;
map { $mean += $_ } @acc;
$mean /= (scalar @acc);
$mean  = sprintf("%.2f", $mean);

print $train_err_fh "Average accuracy: $mean\n";
print $out_fh (join("\t", ($classifier, $mean, $acc[0], $acc[$#acc/4], $acc[$#acc/2], $acc[-$#acc/4], $acc[-1])), "\n");

remove_tree($trainDir, $testDir, $result_file);

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
