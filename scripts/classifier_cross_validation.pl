#!/usr/bin/env perl

=head1 Cross-validate a SciKit classifier/predictor within a "Problem Directory"

    classifier_cross_validation -d probDir -f fraction -c classifier -s sweeps -o output -e log/errors

Mandatory input: A "Problem Directory." (Usual source is build_matrix)

Output: An estimate of the classifier/predictor's accuracy.
Leaves detailed reports in subdirectory ProbDir/Classifiers/ClassifierType/.

=head2 Parameters

=over 4

=item ProbDir

The name of a "Problem Directory" containing the following data:

    ProbDir
        X
        y
        col.h (optional)
        row.h (optional)

X is a tab-seperated matrix of numeric input values.
Each row is an instance; each column is an attribute.

y is a vector of output target values.

row.h and col.h are "human-readable" labels for the matrix/vector indices.
Each line is a tab-seperated triple: (indexNumber, shortLabel, optional_LongLabel).

=back

=cut

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


use ScriptUtils;
my $opt = ScriptUtils::Opts( '',
			     [ 'probDir|d=s',    'Problem Directory',                      {} ],
			     [ 'fraction|f=f',   'Test fraction during cross-validation',  {default => 0.2} ],
			     [ 'sweeps|s=i',     'Number of cross-validation sweeps',      {default =>  10} ],
			     [ 'classifier|c=s', 'Type of classifier or predictor',        {default => 'RandomForestClassifier'} ],
			     [ 'output|o=s',     'Output file (D:STDOUT)',                 {} ],
			     [ 'error|e=s',      'Log/Error file (D:STDERR)',              {} ],
			     [ 'tmpdir|t=s',     'Location for Temporary Directory',       {default => $FIG_Config::temp} ],
			     [ 'verbose|v',      'Verbose mode',                           {} ],
			     [ 'debug|D',        'Debug mode (prevents cleanup)',          {} ],
    );

my $probDir    = $opt->probdir;
my $fraction   = $opt->fraction;
my $num_sweeps = $opt->sweeps;
my $classifier = $opt->classifier;
my $output     = $opt->output;
my $error      = $opt->error;
my $tmpdir     = $opt->tmpdir;
our $verbose   = $opt->verbose || $ENV{VERBOSE};
our $debug     = $opt->debug   || $ENV{DEBUG};

unless ($probDir && $fraction && $classifier) {
    if (@ARGV == 3) {
        ($probDir, $fraction, $classifier) = @ARGV;
    }
    else {
        die "Missing arguments";
    }
}
die "Problem Directory '$probDir' does not exist" unless (-d $probDir);


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
my $trainDir    = $tmpdir . '/tmp.' . basename($probDir) . q(.train.) . $$;
my $testDir     = $tmpdir . '/tmp.' . basename($probDir) . q(.test.) . $$;
my $result_file = $tmpdir . q(/tmp.apply.result.) . $$;


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

    print $err_fh "Doing split\n" if $verbose;
    &run_safe([ "split_probdir", $probDir, $trainDir, $testDir, $fraction ],
              \undef, \undef, $err_fh
        ) || die "Split failed on sweep=$sweep: $?, $!";

    print $err_fh "Doing train\n" if $verbose;
    &run_safe([ "train_classifier", $trainDir, $classifier ],
              \undef, \undef, $train_err_fh
        ) || die "Training failed on sweep=$sweep: $?, $!";

    print $err_fh "Doing apply\n" if $verbose;
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
    my $acc = sprintf("%.2f", $num_instances 
		      ? 100.0 * $num_right / $num_instances
		      : 0.0
	);
    
    push @acc, $acc;
    print $train_err_fh "$acc\n\n";
}
@acc = map { sprintf('%.2f', $_) } sort { $a <=> $b } @acc;


my $mean = 0;
my $num_acc = (scalar @acc);   #...should equal $num_sweeps, but just in case...
map { $mean += $_ } @acc;
$mean /= $num_acc;
$mean  = sprintf("%.2f", $mean);

my ($Qmin, $Q1, $Q2, $Q3, $Qmax) = ($acc[0], $acc[$num_acc/4], $acc[$num_acc/2], $acc[-1-$num_acc/4], $acc[-1]);

my $trimean = (0.25) * ($Q1 + (2.0)*$Q2 + $Q3);   #...Tukey's Trimean (robust mean estimator)...
my $iqr     = ($Q3 - $Q1);                        #...Interquartile Range (robust range estimator)...

print $train_err_fh "Average accuracy: $mean\n";
print $out_fh (join("\t", ($classifier, $mean, $Qmin, $Q1, $Q2, $Q3, $Qmax, $trimean, $iqr)), "\n");

remove_tree($trainDir, $testDir, $result_file) unless $debug;

exit(0);



sub run_safe {
    my ( $args, $in_fh, $out_fh, $err_fh ) = @_;
    print $err_fh Dumper($args, $in_fh, $out_fh, $err_fh) if $debug;

    if (my $rc = run3( $args, $in_fh, $out_fh, $err_fh )) {
    	print $err_fh "//\n\n" if $debug;
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
