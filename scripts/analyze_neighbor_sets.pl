=head1 Create Alignment Files for Neighbor Sets

    analyze_neighbor_sets.pl [options] setDir

This script creates the alignment file for a single neighbor set produced by L<create_neighbor_sets.pl>. The output can then be
used as input to the MEGA-X alignment program, and ultimately to build a tree. The alignment will be designed as multi-sited,
using the universal roles for a particular taxonomic grouping found in the B<CheckG> directory. In addition, L<GenomeTypeObject>
files will be built for all the genomes in the C<GTOs> subdirectory.

=head2 Parameters

The positional parameter is the input directory (into which the output file will be placed). The input directory should contain
a C<rep.genomes> file containing the genome IDs. It will gain a C<prots.meg> file to be used as input to the alignment program.

Additional command-line options are those given below.

=over 4

=item domain

The taxonomic grouping ID to be used to compute the universal roles. The default is C<2>, indicating Bacteria.

=item roleFile

The C<roles.in.subsystems> file containing the role IDs, checksums, and names for the stable roles. The default is the
file of that name in the function predictors directory, or the one in the SEEDtk global directory if the predictors do
not have one.

=item checkDir

The directory containing the EvalG universal role data (see L<EvalCom/Tax>). The default is C<CheckG> in the SEEDtk global
data directory.

=item alignOpts

The name of the MAO file for the MUSCLE alignment in META-X. The default is C<muscle_align_protein.mao> in the SEEDtk global
data directory.

=item buildOpts

The name of the MAO file for the tree build in META-X. The default is C<infer_ME_nucleotide.mao> in the SEEDtk global
data directory.

=item clear

If specified, the GTO directory will be cleared before starting; otherwise, if a GTO already exists, it will be reloaded instead
of being built from PATRIC.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use EvalCom::Tax;
use Stats;
use File::Copy::Recursive;
use GenomeTypeObject;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('setDir',
            ['domain=i', 'taxonomic grouping to use for universal roles', { default => 2 }],
            ['checkDir=s', 'completeness data directory', { default => "$FIG_Config::p3data/CheckG" }],
            ['roleFile|rolefile|R=s', 'role definition file', { default => "$FIG_Config::p3data/roles.in.subsystems" }],
            ['alignOpts|alignopts|ao=s', 'alignment options file', { default => "$FIG_Config::global/muscle_align_protein.mao" }],
            ['buildOpts|buildopts|bo=s', 'tree-build options file', { default => "$FIG_Config::global/infer_ME_nucleotide.mao" }],
            ['clear', 'clear GTO directory before starting']
        );
# Get access to PATRIC.
my $p3 = P3DataAPI->new();
# Get the working directory.
my ($setDir) = @ARGV;
if (! $setDir) {
    die "No input directory specified.";
} elsif (! -s "$setDir/rep.genomes") {
    die "$setDir does not appear to be a neighborhood directory (run create_neighbor_sets first).";
}
# Insure we have the GTO directory.
if (! -d "$setDir/GTOs") {
    print "Creating GTO directory.\n";
    File::Copy::Recursive::pathmk("$setDir/GTOs") || die "Could not create GTO directory for $setDir: $!";
} elsif ($opt->clear) {
    print "Erasing GTO directory.\n";
    File::Copy::Recursive::pathempty("$setDir/GTOs") || die "Could not erase GTO directory for $setDir: $!";
}
my $stats = Stats->new();
# Create the genome checker.
my $checker = EvalCom::Tax->new($opt->checkdir, rolesInSubsystems => $opt->rolefile, logH => \*STDOUT, stats => $stats);
# Get the domain information.
my $domain = $opt->domain;
my ($taxName, $roleHash) = $checker->taxon_data($domain);
if (! $taxName) {
    die "Invalid domain ID $domain.";
} else {
    my $count = scalar(keys %$roleHash);
    print "Using $count universal roles for $domain: $taxName.\n";
}
# Get the other options.
my $alignOpts = $opt->alignopts;
my $buildOpts = $opt->buildopts;
# Read in all the genomes in this neighborhood.
print "Processing $setDir.\n";
open(my $ih, "<$setDir/rep.genomes") || die "Could not open genome list: $!";
my %genomes;
while (! eof $ih) {
    my $line = <$ih>;
    $stats->Add(lineIn => 1);
    if ($line =~ /(\d+\.\d+)/) {
        $genomes{$1} = 1;
        $stats->Add(genomeIn => 1);
    }
}
close $ih;
my $gTotal = scalar(keys %genomes);
print "$gTotal genomes read from $setDir.\n";
# Create the GTOs. Before storing each one, we parse out the universal proteins to eventually
# put into the output file. The following hash will contain all the role proteins, sub-hashed by genome ID.
# The theory is that will be exactly one of each of these roles in each genome, but the truth is slightly
# different. We keep the longest version of each protein.
my %roleProteins;
for my $genome (sort keys %genomes) {
    print "Processing $genome.\n";
    my $gtoFile = "$setDir/GTOs/$genome.gto";
    my $gto;
    my $gtoSaved;
    if (-s $gtoFile) {
        $stats->Add(gtoRead => 1);
        $gto = GenomeTypeObject->create_from_file($gtoFile);
        $gtoSaved = 1;
    } else {
        $gto = $p3->gto_of($genome);
        if (! $gto) {
            die "$genome not found in PATRIC.";
        }
        $stats->Add(gtoFetched => 1);
    }
    # Get all the features.
    my $count = 0;
    my $fidList = $gto->{features};
    for my $fid (@$fidList) {
        $stats->Add(featureChecked => 1);
        # Get the protein.
        my $protSeq = $fid->{protein_translation};
        if (! $protSeq) {
            $stats->Add(featureNoProtein => 1);
        } else {
            my $protLen = length $protSeq;
            my @roles = $checker->roles_of_function($fid->{function});
            for my $role (@roles) {
                $stats->Add(roleFound => 1);
                if ($roleHash->{$role}) {
                    $stats->Add(roleUniversal => 1);
                    # Here we have a universal role.
                    my $oldSeq = $roleProteins{$role}{$genome};
                    if (! $oldSeq || length($oldSeq) < $protLen) {
                        # This is a new protein or it is longer, so save it.
                        $roleProteins{$role}{$genome} = $protSeq;
                        $stats->Add(proteinStored => 1);
                        $count++;
                    }
                }
            }
        }
    }
    print "$count proteins stored.\n";
    # Now we have pulled all the useful proteins from the GTO. Write it to disk if we did not
    # load it from there.
    if (! $gtoSaved) {
        $gto->destroy_to_file("$setDir/GTOs/$genome.gto");
        $stats->Add(gtoOut => 1);
    }
}
# Now we write the .meg file with all the proteins in it. This can be used by MEGA to create the tree.
# It requires we align everything first.
print "Writing FASTA files.\n";
my $protFile = "$setDir/align.meg";
open(my $mh, ">$protFile") || die "Could not open protein output file: $!";
print $mh "#mega\n!Title Protein Alignments for $setDir;\n!Format DataType=Protein indel=-;";
# Loop through the roles.
for my $role (sort keys %roleProteins) {
    my $roleGenomes = $roleProteins{$role};
    # Insure this role is fully populated.
    my $gCount = scalar(keys %$roleGenomes);
    if ($gCount < $gTotal) {
        print "Role $role has only $gCount genomes-- skipped.\n";
        $stats->Add(roleSkipped => 1);
    } else {
        print "Writing $role.fasta\n";
        open(my $oh, ">$setDir/$role.fasta") || die "Could not open $role output: $!";
        $stats->Add(roleOut => 1);
        for my $genome (sort keys %$roleGenomes) {
            print $oh ">$genome\n$roleGenomes->{$genome}\n";
            $stats->Add(proteinOut => 1);
        }
        close $oh;
        # Create the alignment. Note we must delete any old file because MEGA will not overwrite.
        my $megFile = "$setDir/$role.meg";
        if (-f $megFile) {
            unlink $megFile;
        }
        my $rc = system("megacc -a $alignOpts -d $setDir/$role.fasta -o $megFile -n");
        if ($rc) {
            die "megacc failed with code $rc.";
        } else {
            # Reopen the alignment file.
            open(my $ih, "<$megFile") || die "Could not open $megFile: $!";
            # Throw out the header line.
            my $line = <$ih>;
            print $mh "!Gene=$role;\n";
            while (! eof $ih) {
                $line = <$ih>;
                if (! $line || $line=~ /^!/) {
                    $stats->Add(alignLineSkipped => 1);
                } else {
                    print $mh $line;
                    $stats->Add(alignLineWritten => 1);
                }
            }
        }
    }
}
close $mh;
# Make a tree out of the alignment.
print "Building phylogenetic tree.\n";
my $rc = system("megacc -a $buildOpts -d $protFile -o $setDir/tree.nwk");
if ($rc) {
    die "megacc failed with code $rc.";
}
print "All done.\n" . $stats->Show();
