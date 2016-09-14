=head1 Alexa Signature Families 

    AlexaSignatureFamilies [options] set1 set2 resultTable

This script is used to test background jobs for Alexa. It is also a useful skeleton script
for creating new Alexa jobs. Its sole function is to wait for a given number of 2-second
intervals.

=head2 Parameters

The positional parameters are the two genome set names and the output table name.

The typical command-line parameters required by Alexa jobs are supported. These are documented in the L<Job> object.

The following command-line parameters are supported by this script.

=over 4

=item minIn

Minimum fraction of genomes in set 1 that must contain a signature family.

=item maxOut

Maximum fraction of genomes in set 2 that may contain a signature family.

=back

=cut

use strict;
use Job_Config;
use Job;
use P3Signatures;

# Parse the command line and create the job object.
my $jobObject = Job->new('set1 set2 resultTable',
        ['minIn=f', 'minimum fraction of set 1 genomes that must contain a signature family', { default => 0.9 }],
        ['maxOut=f', 'maximum fraction of set 2 genomes that may contain a signature family', { default => 0.1 }]
        );
# This variable is needed outside the eval block.
my $familyCount;
# Protect against errors.
eval {
    # Get the command-line options.
    my $minIn = $jobObject->opt->minin;
    my $maxOut = $jobObject->opt->maxout;
    # Get the set/table names.
    my ($gs1, $gs2, $result) = @ARGV;
    if (! $gs1 || ! $gs2) {
        die "You must specify two input sets.";
    } elsif (! $result) {
        die "You must specify a result table.";
    }
    # Get the sets.
    my $list1 = $jobObject->ReadSet($gs1);
    my $list2 = $jobObject->ReadSet($gs2);
    my $count = (scalar @$list1) + (scalar @$list2);
    $jobObject->Progress("Read $count genomes from sets $gs1 and $gs2.");
    # Process the sets.
    my $dataH = P3Signatures::Progress($list1, $list2, $minIn, $maxOut, $jobObject);
    $familyCount = scalar keys %$dataH;
    $jobObject->Progress("Spooling families to output.");
    $jobObject->StoreTable($result, "signature families in $gs1 but not $gs2", [qw(family.family_id count1 count2 family.product)], $dataH);
};
if ($@) {
    $jobObject->Fail("ERROR computing signatures: $@");
} else {
    $jobObject->Finish("Signatures found: $familyCount.");
}
