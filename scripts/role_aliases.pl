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
use Stats;
use RoleParse;


=head1 Process Aliases in a roles.in.subsystems File

    role_aliases.pl [ options ] aliasFile

This script reads a C<roles.in.subsystems> file and a file of role aliases and produces a new C<roles.unified> file in which the aliases have the same
role ID.  This latter file can be used for completeness processing.

=head2 Parameters

The positional parameter is the name of the role alias file.  This file should be tab-delimited, with each line containing a role and its alias.

The command-line options are the following

=over 4

=item input

The name of the C<roles.in.subsystems> file.  The default is C<roles.in.subsystems> in the SEEDtk global directory.

=item output

The name of the output file.  The default is C<roles.unified> in the SEEDtk global directory.

=back

=cut

$| = 1;
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('aliasFile',
        ['input|i=s', 'input role file', { default => "$FIG_Config::p3data/roles.in.subsystems" }],
        ['output|o=s', 'output role file', { default => "$FIG_Config::p3data/roles.unified" }],
        );
# Create a statistics object.
my $stats = Stats->new();
# First, read the aliases.
my ($aliasFile) = @ARGV;
if (! $aliasFile) {
    die "Missing required alias file name.\n";
} elsif (! -s $aliasFile) {
    die "Alias file $aliasFile not found or invalid.";
}
open(my $ah, '<', $aliasFile) || die "Could not open alias file: $!";
my %aliasMap;
while (! eof $ah) {
    my $line = <$ah>;
    chomp $line;
    # Compute the checksums.
    my ($c1, $c2) = map { RoleParse::Checksum($_) } split /\t/, $line;
    $aliasMap{$c1} = $c2;
    $aliasMap{$c2} = $c1;
    if (scalar(keys %aliasMap) % 2 != 0) {
        die "Funny error at $line.\n";
    }
}
print scalar(keys %aliasMap) . " alias mappings stored.\n";
# Now read and copy the input file.
close $ah;
open(my $ih, '<', $opt->input) || die "Could not open input file: $!";
open(my $oh, '>', $opt->output) || die "Could not open output file: $!";
# This hash maps a checksum to its new ID.
my %idMap;
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($id, $check, $name) = split /\t/, $line;
    if ($idMap{$check}) {
        print "Changing $id to $idMap{$check}\n";
        $id = $idMap{$check};
    } elsif ($aliasMap{$check}) {
        # Here this role has an alias.  Insure our ID is associated with it.
        $idMap{$aliasMap{$check}} = $id;
    }
    print $oh "$id\t$check\t$name\n";
}
print "All done.\n";

