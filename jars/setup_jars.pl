=head1 Set Up Jar Files

This command is run once to set up the jar files in a directory for easy calling.
It expects to be invoked from inside the jar directory, so cd to that directory
and use

    perl setup.pl binDir logDir

where I<binDir> is a directory on your PATH where you have execution privileges
and I<logDir> is where you want the log files to appear.

For each command, the script will create a shell file to invoke the jar with the
appropriate tuning parameters.

=cut

use strict;
use Cwd;
use File::Spec;

# Compute the operating system.
my $winMode = ($^O eq 'MSWin32');
my $delim = ($winMode ? ';' : ':');
my ($binDir, $logDir) = @ARGV;
if (! $binDir) {
    die "No bin directory specified.";
}
my @paths = split $delim, $ENV{PATH};
if (! grep { $_ eq $binDir } @paths) {
    die "Bin directory must be in the PATH.";
}
# Compute the absolute path to the jars.
my $jarBase = File::Spec::rel2abs(cwd());
# Verify we have a log directory.
if (! $logDir) {
    die "No log directory specified.";
} elsif (! -d $logDir) {
    die "Log directory not found or invalid.";
}
if (! opendir(my $dh, $jarBase)) {
    die "Cannot open this directory: $!";
} else {
    # Find the jar files.
    my @jars = grep { $_ =~ /\.jar$/ } readdir $dh;
    closedir $dh;
    for my $jarFile (@jars) {
        # Extract the jar name.
        if ($jarFile =~ /(.+)\.jar/) {
            my $javaName = $1;
            # Check for a parm file.
            my $parms = "";
            if (open(my $ph, '<', "$jarBase/$javaName.txt")) {
                $parms = <$ph>;
                chomp $parms;
                close $ph;
            }
            print "Building java command $javaName.\n";
            my $command = "java $parms -Dlogback.configurationFile=$jarBase/logback.xml -Dlogfile.name=$logDir/$javaName -jar $jarBase/$jarFile";
            # Create the appropriate executable file.
            if ($winMode) {
                open(my $jh, '>', "$binDir/$javaName.cmd") || die "Could not open $javaName command file: $!";
                print $jh "\@echo off\n";
                print $jh "$command \%\*\n";
                close $jh;
            } else {
                open(my $jh, '>', "$binDir/$javaName") || die "Could not open $javaName command file: $!";
                print $jh "$command \$\@\n";
                close $jh;
                chmod 0755, "$binDir/$javaName";
            }
        }
    }
}


