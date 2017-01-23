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

# $#train_roles = 0;
Proc::ParallelLoop::pareach(\@train_roles, \&worker, {Max_Workers => 32});


open($tmp_fh, '<', $tmp_file)
    or die "Could not read-open '$tmp_file'";
my @results = <$tmp_fh>;
close($tmp_fh);

print STDOUT @results;
@results = map { chomp; [ split(/\t/, $_) ] } @results;


my $num_found     = 0;
my $num_correct   = 0;
my $num_incorrect = 0;
my $num_under     = 0;
my $num_over      = 0;
my @mult_over     = ();
foreach my $result (@results) {
    if ($result->[2] != 0.0) {
        ++$num_found;
    }

    if ($result->[1] == $result->[2]) {
        ++$num_correct;
    }
    else {
        ++$num_incorrect;

        #..."1" == "predicted", "2" == "actual"...
        if ($result->[2] > $result->[1]) {
            ++$num_over;
            push @mult_over, $result->[2];
        }
        elsif ($result->[2] < $result->[1]) {
            ++$num_under;
        }
        else {
            warn "This can't happen !!!"
        }
    }
}

print STDERR "\n";
print STDERR ("Pct_roles=\t",     sprintf("%0.1f%%", (100.0 * $num_found   / $num_roles)), "\n");
print STDERR ("Consistency=\t",   sprintf("%0.1f%%", (100.0 * $num_correct / $num_roles)), "\n");
if ($num_incorrect > 0) {
    print STDERR ("Pct_Under=\t", sprintf("%0.1f%%", (100.0 * $num_under   / $num_incorrect)), "\n");
    print STDERR ("Pct_Over=\t",  sprintf("%0.1f%%", (100.0 * $num_over    / $num_incorrect)), "\n");
    print STDERR ("Median_mult_over=\t", $mult_over[$#mult_over/2], "\n");
}
exit(0);



sub worker {
    my ($role) = @_;
#   print STDERR "Processing '$role'\n";

    my $accuracy = &SeedUtils::file_read("$trainDir/Predictors/$role/Classifiers/$classifier/accuracy");
    chomp $accuracy;

    my @fields = split /\t/, $accuracy;
    my $worst  = $fields[2];
    my $median = $fields[4];
    if ($worst < 90.0) {
        print STDERR ("$role\tworst-case accuracy=", sprintf("%.2f", $worst), "-- skipping\n");
        return;
    }

    my $col = $train_role_to_index{$role};
    my $tmpDir = "$FIG_Config::temp/tmp.pareach_eval.$role.$$";

    &run_safe([ 'partition_on_column', $testDir, $tmpDir, $col ],
              \undef, \undef, \*STDERR
        );
    die "After partition_on_column into $tmpDir from $testDir on $col.";
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
