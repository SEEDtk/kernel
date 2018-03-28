#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

use IO::File;
use IPC::Run3;
use File::Path qw(make_path rmtree);
require Proc::ParallelLoop;
use FIG_Config;

use SeedUtils;

my $trouble = 0;
my ($testDir, $trainDir, $classifier) = @ARGV;
$classifier ||= 'RandomForestClassifier';      #...Default to "Random Forest"

my @test_roles;
my @test_colsH = &SeedUtils::file_read("$testDir/col.h");
die "Nothing read from '$testDir/col.h'" unless @test_colsH;
my %test_role_to_index  = map { ($_ =~ /^(\d+)\t(\S+)/) ? ($2 => $1) : () } @test_colsH;
my %test_index_to_role  = map { ($_ =~ /^(\d+)\t(\S+)/) ? ($1 => $2) : () } @test_colsH;
map { $test_roles[$_] = $test_index_to_role{$_} } keys %test_index_to_role;
# die Dumper(\@test_roles, \%test_index_to_role);

my @train_roles;
my @train_colsH = &SeedUtils::file_read("$trainDir/col.h");
die "Nothing read from '$trainDir/col.h'" unless @train_colsH;
my %train_role_to_index = map { ($_ =~ /^(\d+)\t(\S+)/) ? ($2 => $1) : () } @train_colsH;
my %train_index_to_role = map { ($_ =~ /^(\d+)\t(\S+)/) ? ($1 => $2) : () } @train_colsH;
map { $train_roles[$_] = $train_index_to_role{$_} } keys %train_index_to_role;
#die Dumper(\@train_roles, \%train_index_to_role);

if (@test_roles != @train_roles) {
    die "ERROR: col.h size-mismatch between '$testDir' and '$trainDir'";

    for (my $i=0; $i < @train_roles; ++$i) {
        if ($test_roles[$i] != $train_roles[$i]) {
            ++$trouble;
            print STDERR "ERROR: mismatch for role '$i' --- '$test_roles[$i]' vs '$train_roles[$i]'\n";
        }
    }
    die "\nAborting due to role mismatches";
}
my $num_roles = @train_roles;
print STDERR "num_roles=$num_roles\n";


my $tmp_fh;
my $tmp_file = "$FIG_Config::temp/tmp.evaluate.$$";
open($tmp_fh, '>', $tmp_file)
    or die "Could not write-open '$tmp_file'";

# $#train_roles = 10;
Proc::ParallelLoop::pareach(\@train_roles, \&worker, {Max_Workers => 32});


open($tmp_fh, '<', $tmp_file)
    or die "Could not read-open '$tmp_file'";
my @results = <$tmp_fh>;
close($tmp_fh);

print STDOUT @results;
@results = map { chomp; [ split(/\t/, $_) ] } @results;


my $num_present   = 0;
my $num_exact     = 0;
my $num_inexact   = 0;
my $num_correct   = 0;
my $num_incorrect = 0;
my $num_under     = 0;
my $num_over      = 0;
my @mult_over     = ();
foreach my $result (@results) {
    my $predicted = $result->[1];
    my $actual    = $result->[2];

    if ($actual > 0.0) {
        ++$num_present;
    }

    if (($predicted > 0.0) xor ($actual > 0.0)) {
        #...only one value is zero...
        ++$num_incorrect;
    }
    else {
        #...either both zero or both nonzero...
        ++$num_correct;
    }

    if ($predicted == $actual) {
        ++$num_exact;
    }
    else {
        ++$num_inexact;

        #..."1" == "predicted", "2" == "actual"...
        if ($predicted < $actual) {
            ++$num_under;
        }
        elsif ($predicted > $actual) {
            ++$num_over;
            push @mult_over, $predicted;
        }
        else {
            warn "This can't happen !!!"
        }
    }
}

print STDERR "\n";
print STDERR ("Pct_roles=\t",          sprintf("%0.1f%%", (100.0 * $num_present / $num_roles)), "\n");
print STDERR ("Coarse_Consistency=\t", sprintf("%0.1f%%", (100.0 * $num_correct / $num_roles)), "\n");
print STDERR ("Fine_Consistency=\t",   sprintf("%0.1f%%", (100.0 * $num_exact   / $num_roles)), "\n");
if ($num_incorrect > 0) {
    print STDERR ("Pct_Under=\t", sprintf("%0.1f%%", (100.0 * $num_under / $num_inexact)), "\n");
    print STDERR ("Pct_Over=\t",  sprintf("%0.1f%%", (100.0 * $num_over  / $num_inexact)), "\n");
    @mult_over = sort { $a <=> $b } @mult_over;
    print STDERR ("Median_mult_over=\t", $mult_over[$#mult_over/2], "\n");
}

#...Cleanup scratch data...
unlink($tmp_file) unless $ENV{DEBUG};

exit(0);



sub worker {
    my ($role) = @_;
#   print STDERR "Processing '$role'\n";

    my $accuracy = &SeedUtils::file_read("$trainDir/Predictors/$role/Classifiers/$classifier/accuracy");
    my ($worst, $median) = (0, 0);
    if (! $accuracy) {
        print STDERR ("$role\tno accuracy data\t-- warning\n");
    } else {
        chomp $accuracy;

        my @fields = split /\t/, $accuracy;
        $worst  = $fields[2];
        $median = $fields[4];
        if ($worst < 90.0) {
            print STDERR ("$role\tworst-case accuracy=", sprintf("%.2f", $worst), "\t-- warning\n");
            # return;
        }
    }
    my $col = $train_role_to_index{$role};
    my $tmpDir = "$FIG_Config::temp/tmp.pareach_eval.$role.$$";

    &run_safe([ 'partition_on_column', $testDir, $tmpDir, $col ],
              \undef, \undef, \*STDERR
        );
    &run_safe([ "apply_classifier", $tmpDir, "$trainDir/Predictors/$role", $classifier ],
              \undef, \undef, \*STDERR
        ) || die "Application failed: $?, $!";

    my $predicted = &SeedUtils::file_read("$tmpDir/Classifiers/$classifier/y.out");
    chomp $predicted;

    my $actual    = &SeedUtils::file_read("$tmpDir/y");
    chomp $actual;

    print $tmp_fh (join("\t", ($role, $predicted, $actual)), "\n");

    rmtree( $tmpDir );

    if ($predicted ne $actual) {
        if ($median > 98.0) {
            print STDERR (join("\t", ("mismatch:", $role, sprintf("%.2f", $worst), sprintf("%.2f", $median), $predicted, $actual)), "\n");
        }
    }

    return;
}



sub run_safe {
    my ( $args, $in_fh, $out_fh, $err_fh ) = @_;
#   print $err_fh Dumper($args, $in_fh, $out_fh, $err_fh);

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
