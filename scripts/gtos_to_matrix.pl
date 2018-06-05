#!/usr/bin/env perl
=head1 GTO Files to Role Matrix

    gtos_to_matrix.pl gtoList probDir roles.in.subsystems roles.to.use

This program extracts the roles from a L<GenomeTypeObject> and computes the role IDs. It produces a file containing the unmapped roles on the standard error output. Three output files are produced: a matrix 'X' of role counts with genomes as rows and roles as columns, a list of GTO names 'row.h', and a list of role IDs used 'col.h'.

=head2 Parameters

The positional parameters are a file containing a list of GTO input files (in json format), the directory in which to place the output files and the name of a file containing the IDs of roles to keep (produced by L<build_role_tables.pl>).

=cut

use strict;
use warnings;
use Data::Dumper;
use GenomeTypeObject;
use SeedUtils;
use ScriptUtils;
use RoleParse;
use Data::Dump qw(pp);
use File::Copy::Recursive;

my $opt = ScriptUtils::Opts('gto_file probDir roles.in.subsystems roles.to.use',
        ['clear', 'overwrite previous results']);

my ($gto_list, $probDir, $roles_in_subsystems, $roles_to_use) = @ARGV;
my (%roleIDs);

if (!-d $probDir)
{
    mkdir($probDir) or die "Could not create probDir '$probDir'";
}

elsif ($opt->clear)
{
    print "Clearing $probDir.\n";

    if (-d "$probDir/Predictors")
    {
        File::Copy::Recursive::pathrmdir("$probDir/Predictors")
                || die "Could not clear Predictors directory: $!\ ";
    }

    my @files = qw(col.h roles.mapped roles.not_mapped row.h X);

    for my $file (@files)
    {
        if (-f "$probDir/$file")
        {
            unlink "$probDir/$file";
        }
    }
}

else
{
    die "ERROR: probDir '$probDir' already exists";
}

die "Input file '$roles_to_use' does not exist" unless (-s $roles_to_use);

open my $gto_fh, "<", "$gto_list"
    or die "Could not read-open '$gto_list'";

open my $X_fh,     '>', "$probDir/X"
        or die "Could not write-open '$probDir/X'";

open my $row_fh,   '>', "$probDir/row.h"
        or die "Could not write-open '$probDir/row.h'";

open my $col_fh,   '>', "$probDir/col.h"
        or die "Could not write-open '$probDir/col.h'";

open my $unmapped_fh,   '>', "$probDir/roles.not_mapped"
        or die "Could not write-open '$probDir/roles.not_mapped'";

open my $roles_fh,   '<', "$roles_to_use"
        or die "Could not read-open '$roles_to_use'";

my %roleMap = map { my ($roleID, $checksum) = split /\t/;
        ($checksum => $roleID)} &SeedUtils::file_read($roles_in_subsystems);

my %mapRoles = reverse %roleMap;

my $role_count = 0;

while(<$roles_fh>)
{
        chomp(my $ID = $_);
        if ($ID =~ /^(\S+)/) {
            my $roleID = $1;
            my $checksum = $mapRoles{$roleID};
            $roleIDs{$checksum} = $roleID;
            print $col_fh "$role_count\t$roleID\n";
            $role_count = $role_count + 1;
        }
}

my $genome_count = 0;

while(<$gto_fh>)
{
        chomp(my $gto = $_);
        die "Input file '$gto' does not exist" unless (-s $gto);
        my @dirs = split(/\//,$gto);
        my $genome = $dirs[-1];

        foreach (@dirs)
        {
                if (m/([0-9]+\.[0-9]+)/)
                {
                        $genome = $1;
                }
        }

        print $row_fh "$genome_count\t$genome\n";
        $genome_count = $genome_count + 1;

        my $proc = GenomeTypeObject->new({file => $gto});

        my @CDSs = map { [ $_->{id}, $_->{function} ] } grep {
                ($_->{type} eq q(CDS)) || ($_->{id} =~ m{\.peg\.}) } $proc->features();

        my %counts;

        foreach my $checksum (keys %roleIDs)
        {
                $counts{$roleIDs{$checksum}} = 0;
        }

        foreach my $cds (@CDSs)
        {
            my ($fid, $func) = @$cds;
            my @roles = &SeedUtils::roles_of_function($func);

            foreach my $role (@roles)
            {
                # Compute the role's checksum.
                my $checksum = RoleParse::Checksum($role);
                # Compute the ID for this checksum.
                my $id = $roleIDs{$checksum};

                if ($id)
                {
                    $counts{$id}++;
                    if ($id eq '16sRrnaNMeth2') {
                        print "$id found in $fid: $func.\n";
                    }
                }

                else
                {
                    print $unmapped_fh "$gto role not mapped:\t$fid\t$role\n";
                }
            }
        }

        my @role_counts;

        for my $ID (sort keys %counts)
        {
            push(@role_counts, $counts{$ID});
        }

        print $X_fh join("\t", @role_counts) . "\n";
}
