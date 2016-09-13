=head1 Alexa Background Test Job

    AlexaTest [options]

This script is used to test background jobs for Alexa. It is also a useful skeleton script
for creating new Alexa jobs. Its sole function is to wait for a given number of 2-second
intervals.

=head2 Parameters

There are no positional parameters.

The typical command-line parameters required by Alexa jobs are supported. These are documented in the L<Job> object.

The following command-line parameters are supported by this script.

=over 4

=item time

Number of intervals to wait.

=back

=cut

use strict;
use Job_Config;
use Job;

# Parse the command line and create the job object.
my $jobObject = Job->new('', ['time=i', 'number of intervals to wait', { default => 10 }]);
# This variable is needed outside the EVAL block.
my $count = 0;
# Protect against errors.
eval {
    # Get the command-line options.
    my $time = $jobObject->opt->time;
    # Start the loop.
    for (my $i = 1; $i <= $time; $i++) {
        $jobObject->Progress("Starting sleep interval $i.");
        sleep 2;
        $count++;
    }
};
if ($@) {
    $jobObject->Fail("ERROR after $count intervals: $@");
} else {
    $jobObject->Finish("$count intervals completed.");
}
