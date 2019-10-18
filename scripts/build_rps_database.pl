#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


use strict;
use warnings;
use FIG_Config;
use ScriptUtils;
use SeedUtils;
use BlastUtils;
use Cwd;
use File::Spec;
use Shrub;

=head1 Build an RPS BLAST Database

    build_rps_database [ options ] db_prefix out_directory

Build an RPS database from a set of alignments. Each alignment is stored in a directory
whose name matches the function being aligned. The alignment itself is in a FASTA file with the suffix
C<ali.trimmed>. The directory must also contain a single-record file named C<ALIGNMENT_LENGTH> that
contains the length of the alignment. The list of directories can be specified via the standard input
(see L<ScriptUtils/ih_options>) or a container directory name can be specified and all subdirectories
will be scanned.

=head2 Parameters

The positional parameters are the prefix name to give to the database created followed by the
name of the output directory where it will be stored.

If C<inDir> is specified, all of the subdirectories of the specified directory are used as input. Otherwise,
the standard input should contain a tab-delimited file, and the names of the input directories should be
specified in the first column.

The command-line options are those found in L<ScriptUtils/ih_options> (standard input) and
L<Shrub/script_options> (database connection) plus the following.

=over 4

=item inDir

If specified, the name of an input directory. All the subdirectories will be used as input.

=back

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('db_prefix out_directory', ScriptUtils::ih_options(),
        Shrub::script_options(),
        ['inDir|d=s', 'directory of input directories (if specified)']
        );
# Get the positional parameters.
my ($db_prefix, $out_directory) = @ARGV;
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Get the input alignments.
my @sets;
if ($opt->indir) {
    # Here we have an input directory.
    my $inDir = $opt->indir;
    opendir(my $dh, $inDir) || die "Could not open $inDir: $!";
    @sets = map { "$inDir/$_" } grep { ! /^\./ && -d "$inDir/$_" } readdir($dh);
} else {
    # Here we have an input file.
    my $ih = ScriptUtils::IH($opt->input);
    @sets = SeedUtils::read_ids($ih);
}
# Verify the output directory.
if ($out_directory) {
    SeedUtils::verify_dir($out_directory);
} else {
    $out_directory = getcwd();
}

# Output the alignment length output file.
my $ali_len_outfile = "$out_directory/$db_prefix.ali.lens";
open(my $ah, ">", $ali_len_outfile)
    || die "Could not write-open $ali_len_outfile: $!";

# @sets = ( "HeatShocProtGrpe" );
# Build structure \@aligns = [ [ 'alignment-id-1', 'optional title', 'align1.fa' ], ... ];
my @aligns;
foreach my $set (@sets) {
    # Get the function name.
    my (undef, undef, $function) = File::Spec->splitpath($set);
    # Check for an alignment.
    my $aliFasta = "$set/$function.ali.trimmed";
    if (! -s $aliFasta) {
        warn "Could not find $aliFasta --- skipping";
    } else {
        # Get the function description.
        my ($desc) = $shrub->GetFlat('Function', 'Function(id) = ?', [$function], 'description');
        if (! $desc) {
            die "Function $function not found in database.";
        }
        push @aligns, [ $function, $desc, $aliFasta ];
        my $ali_len;
        my $len_file = "$set/ALIGNMENT_LENGTH";
        if (-s $len_file) {
            $ali_len = &SeedUtils::file_head( $len_file, 1 );
            chomp $ali_len;
            die "Invalid Alignment-length metafile '$len_file'"
                unless ($ali_len && ($ali_len > 0));
        }
        else {
            die "Alignment-length metafile '$len_file' does not exist or has zero size";
        }
        print $ah "$function\t$ali_len\n";
    }
}

my $dbfile = "$out_directory/$db_prefix.rps.db";
if (-s $dbfile ) {
    unlink $dbfile;
} else {
    BlastUtils::build_rps_db(\@aligns, $dbfile, { pseudo_master => 1, title => $db_prefix } )
        || die "Could not build RPS database";
}

