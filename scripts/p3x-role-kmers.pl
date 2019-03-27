=head1 Find Signature Kmers for Roles

    p3x-role-kmers.pl [options] outFile

This script will find the signature kmers for specified roles.  It takes as input a FASTA file of protein sequences.  The sequence ID
can be anything.  The sequence comment should be a role name.  It will find signature kmers for each role.

=head2 Parameters

The positional parameter should be the output file name.  The output file will be a L<KmerDb> in JSON format.

The standard input can be overridden using the options in L<P3Utils/ih_options>.  This should be the FASTA file of sequences.

The L<Shrub> database can be specified using L<Shrub/script_options>.

Additional command-line options are as follows.

=over 4

=item dna

If specified, the input sequences are treated as DNA instead of proteins.

=item K

The kmer size.  This is C<8> for proteins and C<14> for DNA.

=item nullFile

If specified, a FASTA file containing sequences that should not be in any group.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use FastA;
use KmerDb;
use Shrub;
use RoleParse;
use Stats;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('outFile', Shrub::script_options(), P3Utils::ih_options(),
        ['dna', 'input sequences are DNA, not proteins'],
        ['kmerSize|kmersize|kmer|K=i', 'size of each kmer'],
        ['nullFile=s', 'FASTA file containing contra-indicative sequences']
        );
my $stats = Stats->new();
# Check the parameters.
my ($outFile) = @ARGV;
if (! $outFile) {
    die "No output file name specified.";
}
my $dna = $opt->dna // 0;
my $K = $opt->kmersize || ($dna ? 14 : 8);
my $nullFile = $opt->nullfile;
# Connect to the Shrub.
print "Connecting to database.\n";
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
print "Opening input file.\n";
my $ihandl = P3Utils::ih($opt);
my $ih = FastA->new($ihandl);
# Create the Kmer database.
my $kmerDB = KmerDb->new(kmerSize => $K, mirror => $dna);
# Here we will cache role IDs, keyed by checksum.
my %roles;
# This is the null group.
my $nullID = 'null';
$kmerDB->AddGroup($nullID, 'null');
# Loop through the sequences.  For each sequence, we compute the role ID and add to that group.
print "Reading input sequences.\n";
my $count = 0;
while ($ih->next()) {
    my $role = $ih->comment;
    $stats->Add(seqIn => 1);
    # Compute the role ID.
    my $checksum = RoleParse::Checksum($role);
    my $roleID = $roles{$checksum};
    if (! $roleID) {
        # Here we have a new role.
        ($roleID) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$checksum], 'id');
        if (! $roleID) {
            print "Not found in Shrub: $role\n";
            $stats->Add(roleNotFound => 1);
            $roleID = $nullID;
        } else {
            $kmerDB->AddGroup($roleID, $role);
            $stats->Add(roleAdded => 1);
        }
    } else {
        $stats->Add(roleFound => 1);
    }
    # Add the sequence to the database.
    $kmerDB->AddSequence($roleID, $ih->left());
    $stats->Add(seqAdded => 1);
    $count++;
    print "$count sequences read.\n" if $count % 1000 == 0;
}
print "$count total sequences processed.\n";
# Now process the null file, if any.
if ($nullFile) {
    print "Reading sequences from $nullFile.\n";
    $ih = FastA->new($nullFile);
    $count = 0;
    while ($ih->next()) {
        $kmerDB->AddSequence($nullID, $ih->left());
        $stats->Add(nullIn => 1);
        $stats->Add(seqAdded => 1);
        $count++;
        print "$count sequences read.\n" if $count % 1000 == 0;
    }
    print "$count null sequences processed.\n";
}
# Delete the null group.
print "Clearing null sequences.\n";
my $deleted = $kmerDB->DeleteGroup($nullID);
print "$deleted kmers removed.\n";
# Compute the discriminators.
print "Computing discriminators.\n";
$kmerDB->ComputeDiscriminators();
print "Saving output file.\n";
$kmerDB->Save($outFile);
print "All done.\n" . $stats->Show();
