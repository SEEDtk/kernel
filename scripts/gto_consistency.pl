#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;
use IPC::Run3;
use File::Path qw(make_path rmtree);

my ($genome, $output_dir, $trainDir, $roles_in_subsystems, $roles_to_use) = @ARGV;

if (!-d $output_dir) {
    make_path($output_dir)
        or die "Could not create output-directory: '$output_dir'";
}
else {
    warn "WARNING -- Output-directory already exists: '$output_dir'";
}


my $mapped_roles_fh = new IO::File;
$mapped_roles_fh->open("> $output_dir/roles.mapped")
    or die "Could not write-open '$output_dir/roles.mapped'";

my $unmapped_roles_fh = new IO::File;
$unmapped_roles_fh->open("> $output_dir/roles.not_mapped")
    or die "Could not write-open '$output_dir/roles.not_mapped'";

&run_safe([ 'gto_to_roles', $genome, $roles_in_subsystems],
          \undef, $mapped_roles_fh, $unmapped_roles_fh
    );
$mapped_roles_fh->close;
$unmapped_roles_fh->close;

die "No roles mapped" unless (-s "$output_dir/roles.mapped");


&run_safe([ 'build_matrix', "$output_dir/roles.mapped", "$output_dir/Evaluate",  $roles_to_use ],
          \undef, \undef, "$output_dir/evaluate.log"
    );


&run_safe([ 'evaluate_assignments', "$output_dir/Evaluate", $trainDir, 'RandomForestClassifier' ],
          \undef, "$output_dir/evaluate.out", "$output_dir/evaluate.log"
    );

exit(0);



sub run_safe {
    my ( $args, $in_fh, $out_fh, $err_fh ) = @_;
    print STDERR (join(q( ), @$args), "\n");
#   print STDERR $err_fh Dumper($in_fh, $out_fh, $err_fh);

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

