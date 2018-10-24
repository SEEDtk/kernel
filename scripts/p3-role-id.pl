=head1 Compute Role IDs for PATRIC Features

    p3-role-id.pl [options]

This script takes as input a file of PATRIC features with a product column and appends the role ID to each.  Roles without an ID in
the Shrub will not be output. If a feature has multiple roles, it will get multiple output lines.

=head2 Parameters

There are no positional parameters.

The standard input can be overridden using the options in L<P3Utils/ih_options>.

The Shrub database can be specified using the options in L<Shrub/script_options>.

Additional command-line options are the following.

=over 4

=item product

The index (1-based) or name of the column containing the functional assignment.  The default is C<product>.

=item nohead

If specified, it is presumed the input file has no headers.

=item batchSize

The number of input records to process in a batch (for performance). The default is C<10>.

=back

=cut

use strict;
use P3Utils;
use Shrub;
use RoleParse;
use Data::Dumper;

# Get the command-line options.
my $opt = P3Utils::script_opts('', Shrub::script_options(), P3Utils::ih_options(),
        ['nohead', 'files have no headers'],
        ['product=s', 'index (1-based) or name of the product column', { default => 'product' }],
        ['batchSize=i', 'number of input records per batch', { default => 10 }]
        );
# Connect to the database.
my $shrub = Shrub->new_for_script($opt);
# Open the input file.
my $ih = P3Utils::ih($opt);
# Read the incoming headers.
my ($outHeaders) = P3Utils::process_headers($ih, $opt, 1);
# Find the product column.
my $keyCol = P3Utils::find_column($opt->product, $outHeaders);
# Form the full header set and write it out.
if (! $opt->nohead) {
    push @$outHeaders, 'role';
    P3Utils::print_cols($outHeaders);
}
# We will cache role IDs in here.
my %checksums;
# Loop through the input.
while (! eof $ih) {
    # Get the input batch.
    my $couplets = P3Utils::get_couplets($ih, $keyCol, $opt);
    # Loop through the couplets, producing the roles.
    for my $couplet (@$couplets) {
        my ($function, $line) = @$couplet;
        if (! ($function =~ /^hypothetical/)) {
            # Look for cached roles. Ask for the rest.
            my @checksums = map { RoleParse::Checksum($_) } SeedUtils::roles_of_function($function);
            my (@roles, @parms);
            for my $checksum (@checksums) {
                if ($checksums{$checksum}) {
                    push @roles, $checksums{$checksum}
                } else {
                    push @parms, $checksum;
                }
            }
            for my $parm (@parms) {
                my ($role) = $shrub->GetFlat('Role', 'Role(checksum) = ?', [$parm], 'id');
                if ($role) {
                    $checksums{$parm} = $role;
                    push @roles, $role;
                }
            }
            for my $role (@roles) {
                P3Utils::print_cols([@$line, $role]);
            }
        }
    }
}