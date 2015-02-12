# -*- perl -*-
########################################################################
# Copyright (c) 2003-2006 University of Chicago and Fellowship
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
########################################################################


package Tracer;

    use strict;
    use FIG_Config;
    use base qw(Exporter);
    use vars qw(@EXPORT @EXPORT_OK);
    @EXPORT = qw(Trace T TSetup QTrace Confess MemTrace Cluck Min Max Assert Open OpenDir TICK EmergencyKey ETracing Constrain Insure ChDir Emergency Warn TraceDump IDHASH);
    @EXPORT_OK = qw(GetFile MergeOptions ParseRecord UnEscape Escape PrintLine PutLine);
    use Carp qw(longmess croak carp confess);
    use CGI;
    use Cwd;
    use Digest::MD5;
    use File::Basename;
    use File::Path;
    use File::stat;
    use LWP::UserAgent;
    use Time::HiRes 'gettimeofday';
    use URI::Escape;
    use Time::Local;
    use POSIX qw(strftime);
    use Fcntl qw(:DEFAULT :flock);
    use Data::Dumper;


    #
    # These are made optional in order to facilitate the SAS release
    # that might need Tracer.
    #
    BEGIN {
        eval {
            require FIG_Config;
        };
        if ($@) {
            $FIG_Config::temp = "/tmp";
        }
    }

=head1 Tracing and Debugging Helpers

=head2 Tracing

This package provides simple tracing for debugging and reporting purposes. To use it simply call the
L</TSetup> or L</ETracing> method to set the options and call L</Trace> to write out trace messages.
L</TSetup> and L</ETracing> both establish a I<trace level> and a list of I<categories>. Similarly,
each trace message has a I<trace level> and I<category> associated with it. Only messages whose trace
level is less than or equal to the setup trace level and whose category is activated will
be written. Thus, a higher trace level on a message indicates that the message
is less likely to be seen, while a higher trace level passed to B<TSetup> means more trace messages will
appear.

=head3 Putting Trace Messages in Your Code

To generate a trace message, use the following syntax.

    Trace($message) if T(errors => 4);

This statement will produce a trace message if the trace level is 4 or more and the C<errors>
category is active. There is a special category C<main> that is always active, so

    Trace($message) if T(main => 4);

will trace if the trace level is 4 or more.

If the category name is the same as the package name, all you need is the number. So, if the
following call is made in the B<Sprout> package, it will appear if the C<Sprout> category is
active and the trace level is 2 or more.

    Trace($message) if T(2);

In scripts, where no package name is available, the category defaults to C<main>.

=head3 Custom Tracing

Many programs have customized tracing configured using the L</TSetup> method. This is no longer
the preferred method, but a knowledge of how custom tracing works can make the more modern
L</Emergency Tracing> easier to understand.

To set up custom tracing, you call the L</TSetup> method. The method takes as input a trace level,
a list of category names, and a destination. The trace level and list of category names are
specified as a space-delimited string. Thus

    TSetup('3 errors Sprout ERDB', 'TEXT');

sets the trace level to 3, activates the C<errors>, C<Sprout>, and C<ERDB> categories, and
specifies that messages should be sent to the standard output.

To turn on tracing for ALL categories, use an asterisk. The call below sets every category to
level 3 and writes the output to the standard error output. This sort of thing might be
useful in a CGI environment.

    TSetup('3 *', 'WARN');

In addition standard error and file output for trace messages, you can specify that the trace messages
be queued. The messages can then be retrieved by calling the L</QTrace> method. This approach
is useful if you are building a web page. Instead of having the trace messages interspersed with
the page output, they can be gathered together and displayed at the end of the page. This makes
it easier to debug page formatting problems.

Finally, you can specify that all trace messages be emitted to a file, or the standard output and
a file at the same time. To trace to a file, specify the filename with an output character in front
of it.

    TSetup('4 SQL', ">$fileName");

To trace to the standard output and a file at the same time, put a C<+> in front of the angle
bracket.

    TSetup('3 *', "+>$fileName");

The flexibility of tracing makes it superior to simple use of directives like C<die> and C<warn>.
Tracer calls can be left in the code with minimal overhead and then turned on only when needed.
Thus, debugging information is available and easily retrieved even when the application is
being used out in the field.

=head3 Trace Levels

There is no hard and fast rule on how to use trace levels. The following is therefore only
a suggestion.

=over 4

=item Error 0

Message indicates an error that may lead to incorrect results or that has stopped the
application entirely.

=item Warning 1

Message indicates something that is unexpected but that probably did not interfere
with program execution.

=item Notice 2

Message indicates the beginning or end of a major task.

=item Information 3

Message indicates a subtask. In the FIG system, a subtask generally relates to a single
genome. This would be a big loop that is not expected to execute more than 500 times or so.

=item Detail 4

Message indicates a low-level loop iteration.

=back

The format of trace messages is important because some utilities analyze trace files.
There are three fields-- the time stamp, the category name, and the text.
The time stamp is between square brackets and the category name between angle brackets.
After the category name there is a colon (C<:>) followed by the message text.
If the square brackets or angle brackets are missing, then the trace management
utilities assume that they are encountering a set of pre-formatted lines.

Note, however, that this formatting is done automatically by the tracing functions. You
only need to know about it if you want to parse a trace file.

=head3 Emergency Tracing

Sometimes, you need a way for tracing to happen automatically without putting parameters
in a form or on the command line. Emergency tracing does this. You invoke emergency tracing
from the debug form, which is accessed from the [[DebugConsole]]. Emergency tracing requires
that you specify a tracing key. For command-line tools, the key is
taken from the C<TRACING> environment variable. For web services, the key is taken from
a cookie. Either way, the key tells the tracing facility who you are, so that you control
the tracing in your environment without stepping on other users.

The key can be anything you want. If you don't have a key, the C<SetPassword> page will
generate one for you.

You can activate and de-activate emergency tracing from the debugging control panel, as
well as display the trace file itself.

To enable emergency tracing in your code, call

    ETracing($cgi)

from a web script and

    ETracing()

from a command-line script.

The web script will look for the tracing key in the cookies, and the command-line
script will look for it in the C<TRACING> environment variable. If you are
using the L</StandardSetup> method or a [[WebApplication]], emergency tracing
will be configured automatically.

=cut

# Declare the configuration variables.

my $Destination = "WARN";   # Description of where to send the trace output.
my $TeeFlag = 0;            # TRUE if output is going to a file and to the
                            # standard output
my %Categories = ( main => 1 );
                            # hash of active category names
my @LevelNames = qw(error warn notice info detail);
my $TraceLevel = 0;         # trace level; a higher trace level produces more
                            # messages
my @Queue = ();             # queued list of trace messages.
my $LastCategory = "main";  # name of the last category interrogated
my $LastLevel = 0;          # level of the last test call
my $SetupCount = 0;         # number of times TSetup called
my $AllTrace = 0;           # TRUE if we are tracing all categories.
my $SavedCGI;               # CGI object passed to ETracing
my $CommandLine;            # Command line passed to StandardSetup
my $Confessions = 0;        # confession count
umask 2;                    # Fix the damn umask so everything is group-writable.

=head2 Tracing Methods

=head3 Setups

    my $count = Tracer::Setups();

Return the number of times L</TSetup> has been called.

This method allows for the creation of conditional tracing setups where, for example, we
may want to set up tracing if nobody else has done it before us.

=cut

sub Setups {
    return $SetupCount;
}

=head3 TSetup

    TSetup($categoryList, $target);

This method is used to specify the trace options. The options are stored as package data
and interrogated by the L</Trace> and L</T> methods.

=over 4

=item categoryList

A string specifying the trace level and the categories to be traced, separated by spaces.
The trace level must come first.

=item target

The destination for the trace output. To send the trace output to a file, specify the file
name preceded by a ">" symbol. If a double symbol is used (">>"), then the data is appended
to the file. Otherwise the file is cleared before tracing begins. Precede the first ">"
symbol with a C<+> to echo output to a file AND to the standard output. In addition to
sending the trace messages to a file, you can specify a special destination. C<HTML> will
cause tracing to the standard output with each line formatted as an HTML paragraph. C<TEXT>
will cause tracing to the standard output as ordinary text. C<ERROR> will cause trace
messages to be sent to the standard error output as ordinary text. C<QUEUE> will cause trace
messages to be stored in a queue for later retrieval by the L</QTrace> method. C<WARN> will
cause trace messages to be emitted as warnings using the B<warn> directive.  C<NONE> will
cause tracing to be suppressed.

=back

=cut

sub TSetup {
    # Get the parameters.
    my ($categoryList, $target) = @_;
    # Parse the category list.
    my @categoryData = split m/\s+/, $categoryList;
    # Extract the trace level.
    $TraceLevel = shift @categoryData;
    # Presume category-based tracing until we learn otherwise.
    $AllTrace = 0;
    # Build the category hash. Note that if we find a "*", we turn on non-category
    # tracing. We must also clear away any pre-existing data.
    %Categories = ( main => 1 );
    for my $category (@categoryData) {
        if ($category eq '*') {
            $AllTrace = 1;
        } else {
            $Categories{lc $category} = 1;
        }
    }
    # Now we need to process the destination information. The most important special
    # case is when we're writing to a file. This is indicated by ">" (overwrite) and
    # ">>" (append). A leading "+" for either indicates that we are also writing to
    # the standard output (tee mode).
    if ($target =~ m/^\+?>>?/) {
        if ($target =~ m/^\+/) {
            $TeeFlag = 1;
            $target = substr($target, 1);
        }
        if ($target =~ m/^>[^>]/) {
            # We need to initialize the file (which clears it).
            open TRACEFILE, $target;
            print TRACEFILE "[" . Now() . "] [notice] [Tracer] Tracing initialized.\n";
            close TRACEFILE;
            # Set to append mode now that the file has been cleared.
            $Destination = ">$target";
        } else {
            $Destination = $target;
        }
    } else {
        $Destination = uc($target);
    }
    # Increment the setup counter.
    $SetupCount++;
}

=head3 SetLevel

    Tracer::SetLevel($newLevel);

Modify the trace level. A higher trace level will cause more messages to appear.

=over 4

=item newLevel

Proposed new trace level.

=back

=cut

sub SetLevel {
    $TraceLevel = $_[0];
}

=head3 ParseDate

    my $time = Tracer::ParseDate($dateString);

Convert a date into a PERL time number. This method expects a date-like string
and parses it into a number. The string must be vaguely date-like or it will
return an undefined value. Our requirement is that a month and day be
present and that three pieces of the date string (time of day, month and day,
year) be separated by likely delimiters, such as spaces, commas, and such-like.

If a time of day is present, it must be in military time with two digits for
everything but the hour.

The year must be exactly four digits.

Additional stuff can be in the string. We presume it's time zones or weekdays or something
equally innocuous. This means, however, that a sufficiently long sentence with date-like
parts in it may be interpreted as a date. Hopefully this will not be a problem.

It should be guaranteed that this method will parse the output of the L</Now> function.

The parameters are as follows.

=over 4

=item dateString

The date string to convert.

=item RETURN

Returns a PERL time, that is, a number of seconds since the epoch, or C<undef> if
the date string is invalid. A valid date string must contain a month and day.

=back

=cut

# Universal month conversion table.
use constant MONTHS => {    Jan =>  0, January   =>  0, '01' =>  0,  '1' =>  0,
                            Feb =>  1, February  =>  1, '02' =>  1,  '2' =>  1,
                            Mar =>  2, March     =>  2, '03' =>  2,  '3' =>  2,
                            Apr =>  3, April     =>  3, '04' =>  3,  '4' =>  3,
                            May =>  4, May       =>  4, '05' =>  4,  '5' =>  4,
                            Jun =>  5, June      =>  5, '06' =>  5,  '6' =>  5,
                            Jul =>  6, July      =>  6, '07' =>  6,  '7' =>  6,
                            Aug =>  7, August    =>  7, '08' =>  7,  '8' =>  7,
                            Sep =>  8, September =>  8, '09' =>  8,  '9' =>  8,
                            Oct =>  9, October  =>   9, '10' =>  9,
                            Nov => 10, November =>  10, '11' => 10,
                            Dec => 11, December =>  11, '12' => 11
                        };

sub ParseDate {
    # Get the parameters.
    my ($dateString) = @_;
    # Declare the return variable.
    my $retVal;
    # Find the month and day of month. There are two ways that can happen. We check for the
    # numeric style first. That way, if the user's done something like "Sun 12/22", then we
    # won't be fooled into thinking the month is Sunday.
    if ($dateString =~ m#\b(\d{1,2})/(\d{1,2})\b# || $dateString =~ m#\b(\w+)\s(\d{1,2})\b#) {
        my ($mon, $mday) = (MONTHS->{$1}, $2);
        # Insist that the month and day are valid.
        if (defined($mon) && $2 >= 1 && $2 <= 31) {
            # Find the time.
            my ($hour, $min, $sec) = (0, 0, 0);
            if ($dateString =~ /\b(\d{1,2}):(\d{2}):(\d{2})\b/) {
                ($hour, $min, $sec) = ($1, $2, $3);
            }
            # Find the year.
            my $year;
            if ($dateString =~ /\b(\d{4})\b/) {
                $year = $1;
            } else {
                # Get the default year, which is this one. Note we must convert it to
                # the four-digit value expected by "timelocal".
                (undef, undef, undef, undef, undef, $year) = localtime();
                $year += 1900;
            }
            $retVal = timelocal($sec, $min, $hour, $mday, $mon, $year);
        }
    }
    # Return the result.
    return $retVal;
}

=head3 LogErrors

    Tracer::LogErrors($fileName);

Route the standard error output to a log file.

=over 4

=item fileName

Name of the file to receive the error output.

=back

=cut

sub LogErrors {
    # Get the file name.
    my ($fileName) = @_;
    # Open the file as the standard error output.
    open STDERR, '>', $fileName;
}

=head3 Trace

    Trace($message);

Write a trace message to the target location specified in L</TSetup>. If there has not been
any prior call to B<TSetup>.

=over 4

=item message

Message to write.

=back

=cut

sub Trace {
    # Get the parameters.
    my ($message) = @_;
    # Strip off any line terminators at the end of the message. We will add
    # new-line stuff ourselves.
    my $stripped = Strip($message);
    # Compute the caller information.
    my ($callPackage, $callFile, $callLine) = caller();
    my $callFileTitle = basename($callFile);
    # Check the caller.
    my $callerInfo = ($callFileTitle ne "Tracer.pm" ? " [$callFileTitle $callLine]" : "");
    # Get the timestamp.
    my $timeStamp = Now();
    # Build the prefix.
    my $level = $LevelNames[$LastLevel] || "($LastLevel)";
    my $prefix = "[$timeStamp] [$level] [$LastCategory]$callerInfo";
    # Format the message.
    my $formatted = "$prefix $stripped";
    # Process according to the destination.
    if ($Destination eq "TEXT") {
        # Write the message to the standard output.
        print "$formatted\n";
    } elsif ($Destination eq "ERROR") {
        # Write the message to the error output. Here, we want our prefix fields.
        print STDERR "$formatted\n";
    } elsif ($Destination eq "WARN") {
        # Emit the message to the standard error output. It is presumed that the
        # error logger will add its own prefix fields, the notable exception being
        # the caller info.
        print STDERR "$callerInfo$stripped\n";
    } elsif ($Destination eq "QUEUE") {
        # Push the message into the queue.
        push @Queue, "$formatted";
    } elsif ($Destination eq "HTML") {
        # Convert the message to HTML.
        my $escapedMessage = CGI::escapeHTML($stripped);
        # The stuff after the first line feed should be pre-formatted.
        my @lines = split m/\s*\n/, $escapedMessage;
        # Get the normal portion.
        my $line1 = shift @lines;
        print "<p>$timeStamp $LastCategory $LastLevel: $line1</p>\n";
        if (@lines) {
            print "<pre>" . join("\n", @lines, "</pre>");
        }
    } elsif ($Destination =~ m/^>>/) {
        # Write the trace message to an output file.
        open(TRACING, $Destination) || confess("Tracing open for \"$Destination\" failed: $!");
        # Lock the file.
        flock TRACING, LOCK_EX;
        print TRACING "$formatted\n";
        close TRACING;
        # If the Tee flag is on, echo it to the standard output.
        if ($TeeFlag) {
            print "$formatted\n";
        }
    }
}

=head3 TraceDump

    TraceDump($title, $object);

Dump an object to the trace log. This method simply calls the C<Dumper>
function, but routes the output to the trace log instead of returning it
as a string. The output is arranged so that it comes out monospaced when
it appears in an HTML trace dump.

=over 4

=item title

Title to give to the object being dumped.

=item object

Reference to a list, hash, or object to dump.

=back

=cut

sub TraceDump {
    # Get the parameters.
    my ($title, $object) = @_;
    # Trace the object.
    Trace("Object dump for $title:\n" . Dumper($object));
}

=head3 T

    my $switch = T($category, $traceLevel);

    or

    my $switch = T($traceLevel);

Return TRUE if the trace level is at or above a specified value and the specified category
is active, else FALSE. If no category is specified, the caller's package name is used.

=over 4

=item category

Category to which the message belongs. If not specified, the caller's package name is
used.

=item traceLevel

Relevant tracing level.

=item RETURN

TRUE if a message at the specified trace level would appear in the trace, else FALSE.

=back

=cut

sub T {
    # Declare the return variable.
    my $retVal = 0;
    # Only proceed if tracing is turned on.
    if ($Destination ne "NONE") {
        # Get the parameters.
        my ($category, $traceLevel) = @_;
        if (!defined $traceLevel) {
            # Here we have no category, so we need to get the calling package.
            # The calling package is normally the first parameter. If it is
            # omitted, the first parameter will be the tracelevel. So, the
            # first thing we do is shift the so-called category into the
            # $traceLevel variable where it belongs.
            $traceLevel = $category;
            my ($package, $fileName, $line) = caller;
            # If there is no calling package, we default to "main".
            if (!$package) {
                $category = "main";
            } else {
                my @cats = split /::/, $package;
                $category = $cats[$#cats];
            }
        }
        # Save the category name and level.
        $LastCategory = $category;
        $LastLevel = $traceLevel;
        # Convert it to lower case before we hash it.
        $category = lc $category;
        # Validate the trace level.
        if (ref $traceLevel) {
            Confess("Bad trace level.");
        } elsif (ref $TraceLevel) {
            Confess("Bad trace config.");
        }
        # Make the check. Note that level 0 shows even if the category is turned off.
        $retVal = ($traceLevel <= $TraceLevel && ($traceLevel == 0 || $AllTrace || exists $Categories{$category}));
    }
    # Return the computed result.
    return $retVal;
}

=head3 QTrace

    my $data = QTrace($format);

Return the queued trace data in the specified format.

=over 4

=item format

C<html> to format the data as an HTML list, C<text> to format it as straight text.

=back

=cut

sub QTrace {
    # Get the parameter.
    my ($format) = @_;
    # Create the return variable.
    my $retVal = "";
    # Only proceed if there is an actual queue.
    if (@Queue) {
        # Process according to the format.
        if ($format =~ m/^HTML$/i) {
            # Convert the queue into an HTML list.
            $retVal = "<ul>\n";
            for my $line (@Queue) {
                my $escapedLine = CGI::escapeHTML($line);
                $retVal .= "<li>$escapedLine</li>\n";
            }
            $retVal .= "</ul>\n";
        } elsif ($format =~ m/^TEXT$/i) {
            # Convert the queue into a list of text lines.
            $retVal = join("\n", @Queue) . "\n";
        }
        # Clear the queue.
        @Queue = ();
    }
    # Return the formatted list.
    return $retVal;
}

=head3 Confess

    Confess($message);

Trace the call stack and abort the program with the specified message. When used with
the OR operator and the L</Assert> method, B<Confess> can function as a debugging assert.
So, for example

    Assert($recNum >= 0) || Confess("Invalid record number $recNum.");

Will abort the program with a stack trace if the value of C<$recNum> is negative.

=over 4

=item message

Message to include in the trace.

=back

=cut

sub Confess {
    # Get the parameters.
    my ($message) = @_;
    # Set up the category and level.
    $LastCategory = "(confess)";
    $LastLevel = 0;
    # Trace the call stack.
    Cluck($message);
    # Increment the confession count.
    $Confessions++;
    # Abort the program.
    croak(">>> $message");
}

=head3 Confessions

    my $count = Tracer::Confessions();

Return the number of calls to L</Confess> by the current task.

=cut

sub Confessions {
    return $Confessions;
}


=head3 SaveCGI

    Tracer::SaveCGI($cgi);

This method saves the CGI object but does not activate emergency tracing.
It is used to allow L</Warn> to work in situations where emergency
tracing is contra-indicated (e.g. the wiki).

=over 4

=item cgi

Active CGI query object.

=back

=cut

sub SaveCGI {
    $SavedCGI = $_[0];
}

=head3 Warn

    Warn($message, @options);

This method traces an important message. If an RSS feed is configured
(via I<FIG_Config::error_feed>) and the tracing destination is C<WARN>,
then the message will be echoed to the feed. In general, a tracing
destination of C<WARN> indicates that the caller is running as a web
service in a production environment; however, this is not a requirement.

To force warnings into the RSS feed even when the tracing destination
is not C<WARN>, simply specify the C<Feed> tracing module. This can be
configured automatically when L</StandardSetup> is used.

The L</Cluck> method calls this one for its final message. Since
L</Confess> calls L</Cluck>, this means that any error which is caught
and confessed will put something in the feed. This insures that someone
will be alerted relatively quickly when a failure occurs.

=over 4

=item message

Message to be traced.

=item options

A list containing zero or more options.

=back

The permissible options are as follows.

=over 4

=item noStack

If specified, then the stack trace is not included in the output.

=back

=cut

sub Warn {
    # Get the parameters.
    my $message = shift @_;
    my %options = map { $_ => 1 } @_;
    # Save $@;
    my $savedError = $@;
    # Trace the message.
    Trace($message);
}


=head3 Assert

    Assert($condition1, $condition2, ... $conditionN);

Return TRUE if all the conditions are true. This method can be used in conjunction with
the OR operator and the L</Confess> method as a debugging assert.
So, for example

    Assert($recNum >= 0) || Confess("Invalid record number $recNum.");

Will abort the program with a stack trace if the value of C<$recNum> is negative.

=cut
sub Assert {
    my $retVal = 1;
    LOOP: for my $condition (@_) {
        if (! $condition) {
            $retVal = 0;
            last LOOP;
        }
    }
    return $retVal;
}

=head3 Cluck

    Cluck($message);

Trace the call stack. Note that for best results, you should qualify the call with a
trace condition. For example,

    Cluck("Starting record parse.") if T(3);

will only trace the stack if the trace level for the package is 3 or more.

=over 4

=item message

Message to include in the trace.

=back

=cut

sub Cluck {
    # Get the parameters.
    my ($message) = @_;
    # Trace what's happening.
    Trace("Stack trace for event: $message");
    # Get the stack trace.
    my @trace = LongMess();
    # Convert the trace to a series of messages.
    for my $line (@trace) {
        # Replace the tab at the beginning with spaces.
        $line =~ s/^\t/    /;
        # Trace the line.
        Trace($line);
    }
    # Issue a warning. This displays the event message and inserts it into the RSS error feed.
    Warn($message);
}

=head3 LongMess

    my @lines = Tracer::LongMess();

Return a stack trace with all tracing methods removed. The return will be in the form of a list
of message strings.

=cut

sub LongMess {
    # Declare the return variable.
    my @retVal = ();
    my $confession = longmess("");
    for my $line (split m/\s*\n/, $confession) {
        unless ($line =~ /Tracer\.pm/) {
            # Here we have a line worth keeping. Push it onto the result list.
            push @retVal, $line;
        }
    }
    # Return the result.
    return @retVal;
}

=head3 ETracing

    ETracing($parameter, %options);

Set up emergency tracing. Emergency tracing is tracing that is turned
on automatically for any program that calls this method. The emergency
tracing parameters are stored in a a file identified by a tracing key.
If this method is called with a CGI object, then the tracing key is
taken from a cookie. If it is called with no parameters, then the tracing
key is taken from an environment variable. If it is called with a string,
the tracing key is that string.

=over 4

=item parameter

A parameter from which the tracing key is computed. If it is a scalar,
that scalar is used as the tracing key. If it is a CGI object, the
tracing key is taken from the C<IP> cookie. If it is omitted, the
tracing key is taken from the C<TRACING> environment variable. If it
is a CGI object and emergency tracing is not on, the C<Trace> and
C<TF> parameters will be used to determine the type of tracing.

=item options

Hash of options. The permissible options are given below.

=over 8

=item destType

Emergency tracing destination type to use if no tracing file is found. The
default is C<WARN>.

=item noParms

If TRUE, then display of the saved CGI parms is suppressed. The default is FALSE.

=item level

The trace level to use if no tracing file is found. The default is C<0>.

=back

=back

=cut

sub ETracing {
    # Get the parameter.
    my ($parameter, %options) = @_;
    # Check for CGI mode.
    if (defined $parameter && ref $parameter eq 'CGI') {
        $SavedCGI = $parameter;
    } else {
        $SavedCGI = undef;
    }
    # Check for the noParms option.
    my $noParms = $options{noParms} || 0;
    # Get the default tracing information.
    my $tracing = $options{level} || 0;
    my $dest = $options{destType} || "WARN";
    # Check for emergency tracing.
    my $tkey = EmergencyKey($parameter);
    my $emergencyFile = EmergencyFileName($tkey);
    if (-e $emergencyFile && (my $stat = stat($emergencyFile))) {
        # We have the file. Read in the data.
        my @tracing = GetFile($emergencyFile);
        # Pull off the time limit.
        my $expire = shift @tracing;
        # Convert it to seconds.
        $expire *= 3600;
        # Check the file data.
        my ($now) = gettimeofday;
        if ($now - $stat->mtime <= $expire) {
            # Emergency tracing is on. Pull off the destination and
            # the trace level;
            $dest = shift @tracing;
            my $level = shift @tracing;
            # Insure Tracer is specified.
            my %moduleHash = map { $_ => 1 } @tracing;
            $moduleHash{Tracer} = 1;
            # Set the trace parameter.
            $tracing = join(" ", $level, sort keys %moduleHash);
        }
    }
    # Convert the destination to a real tracing destination.
    $dest = EmergencyTracingDest($tkey, $dest);
    # Setup the tracing we've determined from all the stuff above.
    TSetup($tracing, $dest);
    # Check to see if we're a web script.
    if (defined $SavedCGI) {
        # Yes we are. Trace the form and environment data if it's not suppressed.
        if (! $noParms) {
            TraceParms($SavedCGI);
        }
        # Check for RAW mode. In raw mode, we print a fake header so that we see everything
        # emitted by the script in its raw form.
        if (T(Raw => 3)) {
            print CGI::header(-type => 'text/plain', -tracing => 'Raw');
        }
    }
}

=head3 EmergencyFileName

    my $fileName = Tracer::EmergencyFileName($tkey);

Return the emergency tracing file name. This is the file that specifies
the tracing information.

=over 4

=item tkey

Tracing key for the current program.

=item RETURN

Returns the name of the file to contain the emergency tracing information.

=back

=cut

sub EmergencyFileName {
    # Get the parameters.
    my ($tkey) = @_;
    # Compute the emergency tracing file name.
    return "$FIG_Config::temp/Emergency$tkey.txt";
}

=head3 EmergencyFileTarget

    my $fileName = Tracer::EmergencyFileTarget($tkey);

Return the emergency tracing target file name. This is the file that receives
the tracing output for file-based tracing.

=over 4

=item tkey

Tracing key for the current program.

=item RETURN

Returns the name of the file to contain the trace output.

=back

=cut

sub EmergencyFileTarget {
    # Get the parameters.
    my ($tkey) = @_;
    # Compute the emergency tracing file name.
    return "$FIG_Config::temp/trace$tkey.log";
}

=head3 EmergencyTracingDest

    my $dest = Tracer::EmergencyTracingDest($tkey, $myDest);

This method converts an emergency tracing destination to a real
tracing destination. The main difference is that if the
destination is C<FILE> or C<APPEND>, we convert it to file
output. If the destination is C<DUAL>, we convert it to file
and standard output.

=over 4

=item tkey

Tracing key for this environment.

=item myDest

Destination from the emergency tracing file.

=item RETURN

Returns a destination that can be passed into L</TSetup>.

=back

=cut

sub EmergencyTracingDest {
    # Get the parameters.
    my ($tkey, $myDest) = @_;
    # Declare the return variable.
    my $retVal = $myDest;
    # Process according to the destination value.
    if ($myDest eq 'FILE') {
        $retVal = ">" . EmergencyFileTarget($tkey);
    } elsif ($myDest eq 'APPEND') {
        $retVal = ">>" . EmergencyFileTarget($tkey);
    } elsif ($myDest eq 'DUAL') {
        $retVal = "+>" . EmergencyFileTarget($tkey);
    } elsif ($myDest eq 'WARN') {
        $retVal = "WARN";
    }
    # Return the result.
    return $retVal;
}

=head3 Emergency

    Emergency($key, $hours, $dest, $level, @modules);

Turn on emergency tracing. This method is normally invoked over the web from
a debugging console, but it can also be called by the C<trace.pl> script.
The caller specifies the duration of the emergency in hours, the desired tracing
destination, the trace level, and a list of the trace modules to activate.
For the length of the duration, when a program in an environment with the
specified tracing key active invokes a Sprout CGI script, tracing will be
turned on automatically. See L</TSetup> for more about tracing setup and
L</ETracing> for more about emergency tracing.

=over 4

=item tkey

The tracing key. This is used to identify the control file and the trace file.

=item hours

Number of hours to keep emergency tracing alive.

=item dest

Tracing destination. If no path information is specified for a file
destination, it is put in the FIG temporary directory.

=item level

Tracing level. A higher level means more trace messages.

=item modules

A list of the tracing modules to activate.

=back

=cut

sub Emergency {
    # Get the parameters.
    my ($tkey, $hours, $dest, $level, @modules) = @_;
    # Create the emergency file.
    my $specFile = EmergencyFileName($tkey);
    my $outHandle = Open(undef, ">$specFile");
    print $outHandle join("\n", $hours, $dest, $level, @modules, "");
}

=head3 EmergencyKey

    my $tkey = EmergencyKey($parameter);

Return the Key to be used for emergency tracing. This could be an IP address,
 a session ID, or a user name, depending on the environment.

=over 4

=item parameter

Parameter defining the method for finding the tracing key. If it is a scalar,
then it is presumed to be the tracing key itself. If it is a CGI object, then
the tracing key is taken from the C<IP> cookie. Otherwise, the tracing key is
taken from the C<TRACING> environment variable.

=item RETURN

Returns the key to be used for labels in emergency tracing.

=back

=cut

sub EmergencyKey {
    # Get the parameters.
    my ($parameter) = @_;
    # Declare the return variable.
    my $retVal;
    # Determine the parameter type.
    if (! defined $parameter) {
        # Here we're supposed to check the environment. If that fails, we
        # get the effective login ID.
        $retVal = $ENV{TRACING} || eval { scalar getpwuid($<) } || 100;
    } else {
        my $ptype = ref $parameter;
        if ($ptype eq 'CGI') {
            # Here we were invoked from a web page. Look for a cookie.
            $retVal = $parameter->cookie('IP');
        } elsif (! $ptype) {
            # Here the key was passed in.
            $retVal = $parameter;
        }
    }
    # If no luck finding a key, use the PID.
    if (! defined $retVal) {
        $retVal = $$;
    }
    # Return the result.
    return $retVal;
}


=head3 TraceParms

    Tracer::TraceParms($cgi);

Trace the CGI parameters at trace level CGI => 3 and the environment variables
at level CGI => 4. A self-referencing URL is traced at level CGI => 2.

=over 4

=item cgi

CGI query object containing the parameters to trace.

=back

=cut

sub TraceParms {
    # Get the parameters.
    my ($cgi) = @_;
    if (T(CGI => 2)) {
        # Here we trace the GET-style URL for the script, but only if it's
        # relatively small.
        my $url = $cgi->url(-relative => 1, -query => 1);
        my $len = length($url);
        if ($len < 500) {
            Trace("[URL] $url");
        } elsif ($len > 2048) {
            Trace("[URL] URL is too long to use with GET ($len characters).");
        } else {
            Trace("[URL] URL length is $len characters.");
        }
    }
    if (T(CGI => 3)) {
        # Here we want to trace the parameter data.
        my @names = $cgi->param;
        for my $parmName (sort @names) {
            # Note we skip the Trace parameters, which are for our use only.
            if ($parmName ne 'Trace' && $parmName ne 'TF') {
                my @values = $cgi->param($parmName);
                Trace("[CGI] $parmName = " . join(", ", @values));
            }
        }
        # Display the request method.
        my $method = $cgi->request_method();
        Trace("Method: $method");
    }
    if (T(CGI => 4)) {
        # Here we want the environment data too.
        for my $envName (sort keys %ENV) {
            Trace("[ENV] $envName = $ENV{$envName}");
        }
    }
}

=head3 TraceImages

    Tracer::TraceImages($htmlString);

Trace information about all of an html document's images. The tracing
will be for type "IMG" at level 3. The image's source string
will be displayed. This is generally either the URL of the image or
raw data for the image itself. If the source is too long, only the first 300
characters will be shown at trace level 3. The entire source will be shown,
however, at trace level 4. This method is not very smart, and might catch
Javascript code, but it is still useful when debugging the arcane
behavior of images in multiple browser environments.

=over 4

=item htmlString

HTML text for an outgoing web page.

=back

=cut

sub TraceImages {
    # Only proceed if we're at the proper trace level.
    if (T(IMG => 3)) {
        # For performance reasons we're manipulating $_[0] instead of retrieving the string
        # into a variable called "$htmlString". This is because we expect html strings to be
        # long, and don't want to copy them any more than we have to.
        Trace(length($_[0]) . " characters in web page.");
        # Loop through the HTML, culling image tags.
        while ($_[0] =~ /<img\s+[^>]+?src="([^"]+)"/sgi) {
            # Extract the source string and determine whether or not it's too long.
            my $srcString = $1;
            my $pos = pos($_[0]) - length($srcString);
            my $excess = length($srcString) - 300;
            # We'll put the display string in here.
            my $srcDisplay = $srcString;
            # If it's a data string, split it at the comma.
            $srcDisplay =~ s/^(data[^,]+,)/$1\n/;
            # If there's no excess or we're at trace level 4, we're done. At level 3 with
            # a long string, however, we only show the first 300 characters.
            if ($excess > 0 && ! T(IMG => 4)) {
                $srcDisplay = substr($srcDisplay,0,300) . "\nplus $excess characters.";
            }
            # Output the trace message.
            Trace("Image tag at position $pos:\n$srcDisplay");
        }
    }
}


=head2 File Utility Methods

=head3 GetFile

    my @fileContents = Tracer::GetFile($fileName);

    or

    my $fileContents = Tracer::GetFile($fileName);

Return the entire contents of a file. In list context, line-ends are removed and
each line is a list element. In scalar context, line-ends are replaced by C<\n>.

=over 4

=item fileName

Name of the file to read.

=item RETURN

In a list context, returns the entire file as a list with the line terminators removed.
In a scalar context, returns the entire file as a string. If an error occurs opening
the file, an empty list will be returned.

=back

=cut

sub GetFile {
    # Get the parameters.
    my ($fileName) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Open the file for input.
    my $handle = Open(undef, "<$fileName");
    # Read the whole file into the return variable, stripping off any terminator
    # characters.
    my $lineCount = 0;
    while (! eof $handle) {
        my $line = <$handle>;
        $lineCount++;
        $line = Strip($line);
        push @retVal, $line;
    }
    # Close it.
    close $handle;
    my $actualLines = @retVal;
    Trace("$actualLines lines read from file $fileName.") if T(File => 2);
    # Return the file's contents in the desired format.
    if (wantarray) {
        return @retVal;
    } else {
        return join "\n", @retVal;
    }
}

=head3 PutFile

    Tracer::PutFile($fileName, \@lines);

Write out a file from a list of lines of text.

=over 4

=item fileName

Name of the output file.

=item lines

Reference to a list of text lines. The lines will be written to the file in order, with trailing
new-line characters. Alternatively, may be a string, in which case the string will be written without
modification.

=back

=cut

sub PutFile {
    # Get the parameters.
    my ($fileName, $lines) = @_;
    # Open the output file.
    my $handle = Open(undef, ">$fileName");
    # Count the lines written.
    if (ref $lines ne 'ARRAY') {
        # Here we have a scalar, so we write it raw.
        print $handle $lines;
        Trace("Scalar put to file $fileName.") if T(File => 3);
    } else {
        # Write the lines one at a time.
        my $count = 0;
        for my $line (@{$lines}) {
            print $handle "$line\n";
            $count++;
        }
        Trace("$count lines put to file $fileName.") if T(File => 3);
    }
    # Close the output file.
    close $handle;
}

=head3 ParseRecord

    my @fields = Tracer::ParseRecord($line);

Parse a tab-delimited data line. The data line is split into field values. Embedded tab
and new-line characters in the data line must be represented as C<\t> and C<\n>, respectively.
These will automatically be converted.

=over 4

=item line

Line of data containing the tab-delimited fields.

=item RETURN

Returns a list of the fields found in the data line.

=back

=cut

sub ParseRecord {
    # Get the parameter.
    my ($line) = @_;
    # Remove the trailing new-line, if any.
    chomp $line;
    # Split the line read into pieces using the tab character.
    my @retVal = split /\t/, $line;
    # Trim and fix the escapes in each piece.
    for my $value (@retVal) {
        # Trim leading whitespace.
        $value =~ s/^\s+//;
        # Trim trailing whitespace.
        $value =~ s/\s+$//;
        # Delete the carriage returns.
        $value =~ s/\r//g;
        # Convert the escapes into their real values.
        $value =~ s/\\t/"\t"/ge;
        $value =~ s/\\n/"\n"/ge;
    }
    # Return the result.
    return @retVal;
}

=head3 Merge

    my @mergedList = Tracer::Merge(@inputList);

Sort a list of strings and remove duplicates.

=over 4

=item inputList

List of scalars to sort and merge.

=item RETURN

Returns a list containing the same elements sorted in ascending order with duplicates
removed.

=back

=cut

sub Merge {
    # Get the input list in sort order.
    my @inputList = sort @_;
    # Only proceed if the list has at least two elements.
    if (@inputList > 1) {
        # Now we want to move through the list splicing out duplicates.
        my $i = 0;
        while ($i < @inputList) {
            # Get the current entry.
            my $thisEntry = $inputList[$i];
            # Find out how many elements duplicate the current entry.
            my $j = $i + 1;
            my $dup1 = $i + 1;
            while ($j < @inputList && $inputList[$j] eq $thisEntry) { $j++; };
            # If the number is nonzero, splice out the duplicates found.
            if ($j > $dup1) {
                splice @inputList, $dup1, $j - $dup1;
            }
            # Now the element at position $dup1 is different from the element before it
            # at position $i. We push $i forward one position and start again.
            $i++;
        }
    }
    # Return the merged list.
    return @inputList;
}

=head3 Open

    my $handle = Open($fileHandle, $fileSpec, $message);

Open a file.

The I<$fileSpec> is essentially the second argument of the PERL C<open>
function. The mode is specified using Unix-like shell information. So, for
example,

    Open(\*LOGFILE, '>>/usr/spool/news/twitlog', "Could not open twit log.");

would open for output appended to the specified file, and

    Open(\*DATASTREAM, "| sort -u >$outputFile", "Could not open $outputFile.");

would open a pipe that sorts the records written and removes duplicates. Note
the use of file handle syntax in the Open call. To use anonymous file handles,
code as follows.

    my $logFile = Open(undef, '>>/usr/spool/news/twitlog', "Could not open twit log.");

The I<$message> parameter is used if the open fails. If it is set to C<0>, then
the open returns TRUE if successful and FALSE if an error occurred. Otherwise, a
failed open will throw an exception and the third parameter will be used to construct
an error message. If the parameter is omitted, a standard message is constructed
using the file spec.

    Could not open "/usr/spool/news/twitlog"

Note that the mode characters are automatically cleaned from the file name.
The actual error message from the file system will be captured and appended to the
message in any case.

    Could not open "/usr/spool/news/twitlog": file not found.

In some versions of PERL the only error message we get is a number, which
corresponds to the C++ C<errno> value.

    Could not open "/usr/spool/news/twitlog": 6.

=over 4

=item fileHandle

File handle. If this parameter is C<undef>, a file handle will be generated
and returned as the value of this method.

=item fileSpec

File name and mode, as per the PERL C<open> function.

=item message (optional)

Error message to use if the open fails. If omitted, a standard error message
will be generated. In either case, the error information from the file system
is appended to the message. To specify a conditional open that does not throw
an error if it fails, use C<0>.

=item RETURN

Returns the name of the file handle assigned to the file, or C<undef> if the
open failed.

=back

=cut

sub Open {
    # Get the parameters.
    my ($fileHandle, $fileSpec, $message) = @_;
    # Attempt to open the file.
    my $rv = open $fileHandle, $fileSpec;
    # If the open failed, generate an error message.
    if (! $rv) {
        # Save the system error message.
        my $sysMessage = $!;
        # See if we need a default message.
        if (!$message) {
            # Clean any obvious mode characters and leading spaces from the
            # filename.
            my ($fileName) = FindNamePart($fileSpec);
            $message = "Could not open \"$fileName\"";
        }
        # Terminate with an error using the supplied message and the
        # error message from the file system.
        Confess("$message: $!");
    }
    # Return the file handle.
    return $fileHandle;
}

=head3 FindNamePart

    my ($fileName, $start, $len) = Tracer::FindNamePart($fileSpec);

Extract the portion of a file specification that contains the file name.

A file specification is the string passed to an C<open> call. It specifies the file
mode and name. In a truly complex situation, it can specify a pipe sequence. This
method assumes that the file name is whatever follows the first angle bracket
sequence.  So, for example, in the following strings the file name is
C</usr/fig/myfile.txt>.

    >>/usr/fig/myfile.txt
    </usr/fig/myfile.txt
    | sort -u > /usr/fig/myfile.txt

If the method cannot find a file name using its normal methods, it will return the
whole incoming string.

=over 4

=item fileSpec

File specification string from which the file name is to be extracted.

=item RETURN

Returns a three-element list. The first element contains the file name portion of
the specified string, or the whole string if a file name cannot be found via normal
methods. The second element contains the start position of the file name portion and
the third element contains the length.

=back

=cut
#: Return Type $;
sub FindNamePart {
    # Get the parameters.
    my ($fileSpec) = @_;
    # Default to the whole input string.
    my ($retVal, $pos, $len) = ($fileSpec, 0, length $fileSpec);
    # Parse out the file name if we can.
    if ($fileSpec =~ m/(<|>>?)(.+?)(\s*)$/) {
        $retVal = $2;
        $len = length $retVal;
        $pos = (length $fileSpec) - (length $3) - $len;
    }
    # Return the result.
    return ($retVal, $pos, $len);
}

=head3 OpenDir

    my @files = OpenDir($dirName, $filtered, $flag);

Open a directory and return all the file names. This function essentially performs
the functions of an C<opendir> and C<readdir>. If the I<$filtered> parameter is
set to TRUE, all filenames beginning with a period (C<.>), dollar sign (C<$>),
or pound sign (C<#>) and all filenames ending with a tilde C<~>) will be
filtered out of the return list. If the directory does not open and I<$flag> is not
set, an exception is thrown. So, for example,

    my @files = OpenDir("/Volumes/fig/contigs", 1);

is effectively the same as

    opendir(TMP, "/Volumes/fig/contigs") || Confess("Could not open /Volumes/fig/contigs.");
    my @files = grep { $_ !~ /^[\.\$\#]/ && $_ !~ /~$/ } readdir(TMP);

Similarly, the following code

    my @files = grep { $_ =~ /^\d/ } OpenDir("/Volumes/fig/orgs", 0, 1);

Returns the names of all files in C</Volumes/fig/orgs> that begin with digits and
automatically returns an empty list if the directory fails to open.

=over 4

=item dirName

Name of the directory to open.

=item filtered

TRUE if files whose names begin with a period (C<.>) should be automatically removed
from the list, else FALSE.

=item flag

TRUE if a failure to open is okay, else FALSE

=back

=cut
#: Return Type @;
sub OpenDir {
    # Get the parameters.
    my ($dirName, $filtered, $flag) = @_;
    # Declare the return variable.
    my @retVal = ();
    # Open the directory.
    if (opendir(my $dirHandle, $dirName)) {
        # The directory opened successfully. Get the appropriate list according to the
        # strictures of the filter parameter.
        if ($filtered) {
            @retVal = grep { $_ !~ /^[\.\$\#]/ && $_ !~ /~$/ } readdir $dirHandle;
        } else {
            @retVal = readdir $dirHandle;
        }
        closedir $dirHandle;
    } elsif (! $flag) {
        # Here the directory would not open and it's considered an error.
        Confess("Could not open directory $dirName.");
    }
    # Return the result.
    return @retVal;
}


=head3 Insure

    Insure($dirName, $chmod);

Insure a directory is present.

=over 4

=item dirName

Name of the directory to check. If it does not exist, it will be created.

=item chmod (optional)

Security privileges to be given to the directory if it is created.

=back

=cut

sub Insure {
    my ($dirName, $chmod) = @_;
    if (! -d $dirName) {
        Trace("Creating $dirName directory.") if T(2);
        eval {
            mkpath $dirName;
            # If we have permissions specified, set them here.
            if (defined($chmod)) {
                chmod $chmod, $dirName;
            }
        };
        if ($@) {
            Confess("Error creating $dirName: $@");
        }
    }
}

=head3 ChDir

    ChDir($dirName);

Change to the specified directory.

=over 4

=item dirName

Name of the directory to which we want to change.

=back

=cut

sub ChDir {
    my ($dirName) = @_;
    if (! -d $dirName) {
        Confess("Cannot change to directory $dirName: no such directory.");
    } else {
        Trace("Changing to directory $dirName.") if T(File => 4);
        my $okFlag = chdir $dirName;
        if (! $okFlag) {
            Confess("Error switching to directory $dirName.");
        }
    }
}

=head3 SetPermissions

    Tracer::SetPermissions($dirName, $group, $mask, %otherMasks);

Set the permissions for a directory and all the files and folders inside it.
In addition, the group ownership will be changed to the specified value.

This method is more vulnerable than most to permission and compatability
problems, so it does internal error recovery.

=over 4

=item dirName

Name of the directory to process.

=item group

Name of the group to be assigned.

=item mask

Permission mask. Bits that are C<1> in this mask will be ORed into the
permission bits of any file or directory that does not already have them
set to 1.

=item otherMasks

Map of search patterns to permission masks. If a directory name matches
one of the patterns, that directory and all its members and subdirectories
will be assigned the new pattern. For example, the following would
assign 0664 to most files, but would use 0777 for directories named C<tmp>.

    Tracer::SetPermissions($dirName, 'fig', 01664, '^tmp$' => 01777);

The list is ordered, so the following would use 0777 for C<tmp1> and
0666 for C<tmp>, C<tmp2>, or C<tmp3>.

    Tracer::SetPermissions($dirName, 'fig', 01664, '^tmp1' => 0777,
                                                   '^tmp' => 0666);

Note that the pattern matches are all case-insensitive, and only directory
names are matched, not file names.

=back

=cut

sub SetPermissions {
    # Get the parameters.
    my ($dirName, $group, $mask, @otherMasks) = @_;
    # Set up for error recovery.
    eval {
        # Switch to the specified directory.
        ChDir($dirName);
        # Get the group ID.
        my $gid = getgrnam($group);
        # Get the mask for tracing.
        my $traceMask = sprintf("%04o", $mask) . "($mask)";
        Trace("Fixing permissions for directory $dirName using group $group($gid) and mask $traceMask.") if T(File => 2);
        my $fixCount = 0;
        my $lookCount = 0;
        # @dirs will be a stack of directories to be processed.
        my @dirs = (getcwd());
        while (scalar(@dirs) > 0) {
            # Get the current directory.
            my $dir = pop @dirs;
            # Check for a match to one of the specified directory names. To do
            # that, we need to pull the individual part of the name off of the
            # whole path.
            my $simpleName = $dir;
            if ($dir =~ m!/([^/]+)$!) {
                $simpleName = $1;
            }
            Trace("Simple directory name for $dir is $simpleName.") if T(File => 4);
            # Search for a match.
            my $match = 0;
            my $i;
            for ($i = 0; $i < $#otherMasks && ! $match; $i += 2) {
                my $pattern = $otherMasks[$i];
                if ($simpleName =~ /$pattern/i) {
                    $match = 1;
                }
            }
            # Find out if we have a match. Note we use $i-1 because the loop added 2
            # before terminating due to the match.
            if ($match && $otherMasks[$i-1] != $mask) {
                # This directory matches one of the incoming patterns, and it's
                # a different mask, so we process it recursively with that mask.
                SetPermissions($dir, $group, $otherMasks[$i-1], @otherMasks);
            } else {
                # Here we can process normally. Get all of the non-hidden members.
                my @submems = OpenDir($dir, 1);
                for my $submem (@submems) {
                    # Get the full name.
                    my $thisMem = "$dir/$submem";
                    Trace("Checking member $thisMem.") if T(4);
                    $lookCount++;
                    if ($lookCount % 1000 == 0) {
                        Trace("$lookCount members examined. Current is $thisMem. Mask is $traceMask") if T(File => 3);
                    }
                    # Fix the group.
                    chown -1, $gid, $thisMem;
                    # Insure this member is not a symlink.
                    if (! -l $thisMem) {
                        # Get its info.
                        my $fileInfo = stat $thisMem;
                        # Only proceed if we got the info. Otherwise, it's a hard link
                        # and we want to skip it anyway.
                        if ($fileInfo) {
                            my $fileMode = $fileInfo->mode;
                            if (($fileMode & $mask) != $mask) {
                                # Fix this member.
                                $fileMode |= $mask;
                                chmod $fileMode, $thisMem;
                                $fixCount++;
                            }
                            # If it's a subdirectory, stack it.
                            if (-d $thisMem) {
                                push @dirs, $thisMem;
                            }
                        }
                    }
                }
            }
        }
        Trace("$lookCount files and directories processed, $fixCount fixed.") if T(File => 2);
    };
    # Check for an error.
    if ($@) {
        Confess("SetPermissions error: $@");
    }
}

=head3 GetLine

    my @data = Tracer::GetLine($handle);

Read a line of data from a tab-delimited file.

=over 4

=item handle

Open file handle from which to read.

=item RETURN

Returns a list of the fields in the record read. The fields are presumed to be
tab-delimited. If we are at the end of the file, then an empty list will be
returned. If an empty line is read, a single list item consisting of a null
string will be returned.

=back

=cut

sub GetLine {
    # Get the parameters.
    my ($handle) = @_;
    # Declare the return variable.
    my @retVal = ();
    Trace("File position is " . tell($handle) . ". EOF flag is " . eof($handle) . ".") if T(File => 4);
    # Read from the file.
    my $line = <$handle>;
    # Only proceed if we found something.
    if (defined $line) {
        # Remove the new-line. We are a bit over-cautious here because the file may be coming in via an
        # upload control and have a nonstandard EOL combination.
        $line =~ s/(\r|\n)+$//;
        # Here we do some fancy tracing to help in debugging complicated EOL marks.
        if (T(File => 4)) {
            my $escapedLine = $line;
            $escapedLine =~ s/\n/\\n/g;
            $escapedLine =~ s/\r/\\r/g;
            $escapedLine =~ s/\t/\\t/g;
            Trace("Line read: -->$escapedLine<--");
        }
        # If the line is empty, return a single empty string; otherwise, parse
        # it into fields.
        if ($line eq "") {
            push @retVal, "";
        } else {
            push @retVal, split /\t/,$line;
        }
    } else {
        # Trace the reason the read failed.
        Trace("End of file: $!") if T(File => 3);
    }
    # Return the result.
    return @retVal;
}

=head3 PutLine

    Tracer::PutLine($handle, \@fields, $eol);

Write a line of data to a tab-delimited file. The specified field values will be
output in tab-separated form, with a trailing new-line.

=over 4

=item handle

Output file handle.

=item fields

List of field values.

=item eol (optional)

End-of-line character (default is "\n").

=back

=cut

sub PutLine {
    # Get the parameters.
    my ($handle, $fields, $eol) = @_;
    # Write the data.
    print $handle join("\t", @{$fields}) . ($eol || "\n");
}


=head3 PrintLine

    Tracer::PrintLine($line);

Print a line of text with a trailing new-line.

=over 4

=item line

Line of text to print.

=back

=cut

sub PrintLine {
    # Get the parameters.
    my ($line) = @_;
    # Print the line.
    print "$line\n";
}


=head2 Other Useful Methods

=head3 MergeOptions

    Tracer::MergeOptions(\%table, \%defaults);

Merge default values into a hash table. This method looks at the key-value pairs in the
second (default) hash, and if a matching key is not found in the first hash, the default
pair is copied in. The process is similar to L</GetOptions>, but there is no error-
checking and no return value.

=over 4

=item table

Hash table to be updated with the default values.

=item defaults

Default values to be merged into the first hash table if they are not already present.

=back

=cut

sub MergeOptions {
    # Get the parameters.
    my ($table, $defaults) = @_;
    # Loop through the defaults.
    while (my ($key, $value) = each %{$defaults}) {
        if (!exists $table->{$key}) {
            $table->{$key} = $value;
        }
    }
}

=head3 IDHASH

    my $hash = SHTargetSearch::IDHASH(@keys);

This is a dinky little method that converts a list of values to a reference
to hash of values to labels. The values and labels are the same.

=cut

sub IDHASH {
    my %retVal = map { $_ => $_ } @_;
    return \%retVal;
}

=head3 Pluralize

    my $plural = Tracer::Pluralize($word);

This is a very simple pluralization utility. It adds an C<s> at the end
of the input word unless it already ends in an C<s>, in which case it
adds C<es>.

=over 4

=item word

Singular word to pluralize.

=item RETURN

Returns the probable plural form of the word.

=back

=cut

sub Pluralize {
    # Get the parameters.
    my ($word) = @_;
    # Declare the return variable.
    my $retVal;
    if ($word =~ /s$/) {
        $retVal = $word . 'es';
    } else {
        $retVal = $word . 's';
    }
    # Return the result.
    return $retVal;
}

=head3 Numeric

    my $okFlag = Tracer::Numeric($string);

Return the value of the specified string if it is numeric, or an undefined value
if it is not numeric.

=over 4

=item string

String to check.

=item RETURN

Returns the numeric value of the string if successful, or C<undef> if the string
is not numeric.

=back

=cut

sub Numeric {
    # Get the parameters.
    my ($string) = @_;
    # We'll put the value in here if we succeed.
    my $retVal;
    # Get a working copy of the string.
    my $copy = $string;
    # Trim leading and trailing spaces.
    $copy =~ s/^\s+//;
    $copy =~ s/\s+$//;
    # Check the result.
    if ($copy =~ /^[+-]?\d+$/) {
        $retVal = $copy;
    } elsif ($copy =~ /^([+-]\d+|\d*)[eE][+-]?\d+$/) {
        $retVal = $copy;
    } elsif ($copy =~ /^([+-]\d+|\d*)\.\d*([eE][+-]?\d+)?$/) {
        $retVal = $copy;
    }
    # Return the result.
    return $retVal;
}


=head3 ParseParm

    my $listValue = Tracer::ParseParm($string);

Convert a parameter into a list reference. If the parameter is undefined,
an undefined value will be returned. Otherwise, it will be parsed as a
comma-separated list of values.

=over 4

=item string

Incoming string.

=item RETURN

Returns a reference to a list of values, or C<undef> if the incoming value
was undefined.

=back

=cut

sub ParseParm {
    # Get the parameters.
    my ($string) = @_;
    # Declare the return variable.
    my $retVal;
    # Check for data.
    if (defined $string) {
        # We have some, so split it into a list.
        $retVal = [ split /\s*,\s*/, $string];
    }
    # Return the result.
    return $retVal;
}

=head3 Now

    my $string = Tracer::Now();

Return a displayable time stamp containing the local time. Whatever format this
method produces must be parseable by L</ParseDate>.

=cut

sub Now {
    return DisplayTime(time);
}

=head3 DisplayTime

    my $string = Tracer::DisplayTime($time);

Convert a time value to a displayable time stamp. Whatever format this
method produces must be parseable by L</ParseDate>.

=over 4

=item time

Time to display, in seconds since the epoch, or C<undef> if the time is unknown.

=item RETURN

Returns a displayable time, or C<(n/a)> if the incoming time is undefined.

=back

=cut

sub DisplayTime {
    my ($time) = @_;
    my $retVal = "(n/a)";
    if (defined $time) {
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
        $retVal = _p2($mon+1) . "/" . _p2($mday) . "/" . ($year + 1900) . " " .
                  _p2($hour) . ":" . _p2($min) . ":" . _p2($sec);
    }
    return $retVal;
}

# Pad a number to 2 digits.
sub _p2 {
    my ($value) = @_;
    $value = "0$value" if ($value < 10);
    return $value;
}

=head3 Escape

    my $codedString = Tracer::Escape($realString);

Escape a string for use in a command. Tabs will be replaced by C<\t>, new-lines
replaced by C<\n>, carriage returns will be deleted, and backslashes will be doubled. The
result is to reverse the effect of L</UnEscape>.

=over 4

=item realString

String to escape.

=item RETURN

Escaped equivalent of the real string.

=back

=cut

sub Escape {
    # Get the parameter.
    my ($realString) = @_;
    # Initialize the return variable.
    my $retVal = "";
    # Loop through the parameter string, looking for sequences to escape.
    while (length $realString > 0) {
        # Look for the first sequence to escape.
        if ($realString =~ /^(.*?)([\n\t\r\\])/) {
            # Here we found it. The text preceding the sequence is in $1. The sequence
            # itself is in $2. First, move the clear text to the return variable.
            $retVal .= $1;
            # Strip the processed section off the real string.
            $realString = substr $realString, (length $2) + (length $1);
            # Get the matched character.
            my $char = $2;
            # If we have a CR, we are done.
            if ($char ne "\r") {
                # It's not a CR, so encode the escape sequence.
                $char =~ tr/\t\n/tn/;
                $retVal .= "\\" . $char;
            }
        } else {
            # Here there are no more escape sequences. The rest of the string is
            # transferred unmodified.
            $retVal .= $realString;
            $realString = "";
        }
    }
    # Return the result.
    return $retVal;
}

=head3 UnEscape

    my $realString = Tracer::UnEscape($codedString);

Replace escape sequences with their actual equivalents. C<\t> will be replaced by
a tab, C<\n> by a new-line character, and C<\\> by a backslash. C<\r> codes will
be deleted.

=over 4

=item codedString

String to un-escape.

=item RETURN

Returns a copy of the original string with the escape sequences converted to their actual
values.

=back

=cut

sub UnEscape {
    # Get the parameter.
    my ($codedString) = @_;
    # Initialize the return variable.
    my $retVal = "";
    # Only proceed if the incoming string is nonempty.
    if (defined $codedString) {
        # Loop through the parameter string, looking for escape sequences. We can't do
        # translating because it causes problems with the escaped slash. ("\\t" becomes
        # "\<tab>" no matter what we do.)
        while (length $codedString > 0) {
            # Look for the first escape sequence.
            if ($codedString =~ /^(.*?)\\(\\|n|t|r)/) {
                # Here we found it. The text preceding the sequence is in $1. The sequence
                # itself is in $2. First, move the clear text to the return variable.
                $retVal .= $1;
                $codedString = substr $codedString, (2 + length $1);
                # Get the escape value.
                my $char = $2;
                # If we have a "\r", we are done.
                if ($char ne 'r') {
                    # Here it's not an 'r', so we convert it.
                    $char =~ tr/\\tn/\\\t\n/;
                    $retVal .= $char;
                }
            } else {
                # Here there are no more escape sequences. The rest of the string is
                # transferred unmodified.
                $retVal .= $codedString;
                $codedString = "";
            }
        }
    }
    # Return the result.
    return $retVal;
}

=head3 Percent

    my $percent = Tracer::Percent($number, $base);

Returns the percent of the base represented by the given number. If the base
is zero, returns zero.

=over 4

=item number

Percent numerator.

=item base

Percent base.

=item RETURN

Returns the percentage of the base represented by the numerator.

=back

=cut

sub Percent {
    # Get the parameters.
    my ($number, $base) = @_;
    # Declare the return variable.
    my $retVal = 0;
    # Compute the percent.
    if ($base != 0) {
        $retVal = $number * 100 / $base;
    }
    # Return the result.
    return $retVal;
}

=head3 In

    my $flag = Tracer::In($value, $min, $max);

Return TRUE if the value is between the minimum and the maximum, else FALSE.

=cut

sub In {
    return ($_[0] <= $_[2] && $_[0] >= $_[1]);
}


=head3 Constrain

    my $constrained = Constrain($value, $min, $max);

Modify a numeric value to bring it to a point in between a maximum and a minimum.

=over 4

=item value

Value to constrain.

=item min (optional)

Minimum permissible value. If this parameter is undefined, no minimum constraint will be applied.

=item max (optional)

Maximum permissible value. If this parameter is undefined, no maximum constraint will be applied.

=item RETURN

Returns the incoming value, constrained according to the other parameters.

=back

=cut

sub Constrain {
    # Get the parameters.
    my ($value, $min, $max) = @_;
    # Declare the return variable.
    my $retVal = $value;
    # Apply the minimum constraint.
    if (defined $min && $retVal < $min) {
        $retVal = $min;
    }
    # Apply the maximum constraint.
    if (defined $max && $retVal > $max) {
        $retVal = $max;
    }
    # Return the result.
    return $retVal;
}

=head3 Min

    my $min = Min($value1, $value2, ... $valueN);

Return the minimum argument. The arguments are treated as numbers.

=over 4

=item $value1, $value2, ... $valueN

List of numbers to compare.

=item RETURN

Returns the lowest number in the list.

=back

=cut

sub Min {
    # Get the parameters. Note that we prime the return value with the first parameter.
    my ($retVal, @values) = @_;
    # Loop through the remaining parameters, looking for the lowest.
    for my $value (@values) {
        if ($value < $retVal) {
            $retVal = $value;
        }
    }
    # Return the minimum found.
    return $retVal;
}

=head3 Max

    my $max = Max($value1, $value2, ... $valueN);

Return the maximum argument. The arguments are treated as numbers.

=over 4

=item $value1, $value2, ... $valueN

List of numbers to compare.

=item RETURN

Returns the highest number in the list.

=back

=cut

sub Max {
    # Get the parameters. Note that we prime the return value with the first parameter.
    my ($retVal, @values) = @_;
    # Loop through the remaining parameters, looking for the highest.
    for my $value (@values) {
        if ($value > $retVal) {
            $retVal = $value;
        }
    }
    # Return the maximum found.
    return $retVal;
}

=head3 Strip

    my $string = Tracer::Strip($line);

Strip all line terminators off a string. This is necessary when dealing with files
that may have been transferred back and forth several times among different
operating environments.

=over 4

=item line

Line of text to be stripped.

=item RETURN

The same line of text with all the line-ending characters chopped from the end.

=back

=cut

sub Strip {
    # Get a copy of the parameter string.
    my ($string) = @_;
    my $retVal = (defined $string ? $string : "");
    # Strip the line terminator characters.
    $retVal =~ s/(\r|\n)+$//g;
    # Return the result.
    return $retVal;
}

=head3 Trim

    my $string = Tracer::Trim($line);

Trim all spaces from the beginning and ending of a string.

=over 4

=item line

Line of text to be trimmed.

=item RETURN

The same line of text with all whitespace chopped off either end.

=back

=cut

sub Trim {
    # Get a copy of the parameter string.
    my ($string) = @_;
    my $retVal = (defined $string ? $string : "");
    # Strip the front spaces.
    $retVal =~ s/^\s+//;
    # Strip the back spaces.
    $retVal =~ s/\s+$//;
    # Return the result.
    return $retVal;
}

=head3 Pad

    my $paddedString = Tracer::Pad($string, $len, $left, $padChar);

Pad a string to a specified length. The pad character will be a
space, and the padding will be on the right side unless specified
in the third parameter.

=over 4

=item string

String to be padded.

=item len

Desired length of the padded string.

=item left (optional)

TRUE if the string is to be left-padded; otherwise it will be padded on the right.

=item padChar (optional)

Character to use for padding. The default is a space.

=item RETURN

Returns a copy of the original string with the pad character added to the
specified end so that it achieves the desired length.

=back

=cut

sub Pad {
    # Get the parameters.
    my ($string, $len, $left, $padChar) = @_;
    # Compute the padding character.
    if (! defined $padChar) {
        $padChar = " ";
    }
    # Compute the number of spaces needed.
    my $needed = $len - length $string;
    # Copy the string into the return variable.
    my $retVal = $string;
    # Only proceed if padding is needed.
    if ($needed > 0) {
        # Create the pad string.
        my $pad = $padChar x $needed;
        # Affix it to the return value.
        if ($left) {
            $retVal = $pad . $retVal;
        } else {
            $retVal .= $pad;
        }
    }
    # Return the result.
    return $retVal;
}

=head3 Quoted

    my $string = Tracer::Quoted($var);

Convert the specified value to a string and enclose it in single quotes.
If it's undefined, the string C<undef> in angle brackets will be used
instead.

=over 4

=item var

Value to quote.

=item RETURN

Returns a string enclosed in quotes, or an indication the value is undefined.

=back

=cut

sub Quoted {
    # Get the parameters.
    my ($var) = @_;
    # Declare the return variable.
    my $retVal;
    # Are we undefined?
    if (! defined $var) {
        $retVal = "<undef>";
    } else {
        # No, so convert to a string and enclose in quotes.
        $retVal = $var;
        $retVal =~ s/'/\\'/;
        $retVal = "'$retVal'";
    }
    # Return the result.
    return $retVal;
}

=head3 EOF

This is a constant that is lexically greater than any useful string.

=cut

sub EOF {
    return "\xFF\xFF\xFF\xFF\xFF";
}

=head3 TICK

    my @results = TICK($commandString);

Perform a back-tick operation on a command. If this is a Windows environment, any leading
dot-slash (C<./> will be removed. So, for example, if you were doing

    `./protein.cgi`

from inside a CGI script, it would work fine in Unix, but would issue an error message
in Windows complaining that C<'.'> is not a valid command. If instead you code

    TICK("./protein.cgi")

it will work correctly in both environments.

=over 4

=item commandString

The command string to pass to the system.

=item RETURN

Returns the standard output from the specified command, as a list.

=back

=cut
#: Return Type @;
sub TICK {
    # Get the parameters.
    my ($commandString) = @_;
    # Chop off the dot-slash if this is Windows.
    if ($FIG_Config::win_mode) {
        $commandString =~ s!^\./!!;
    }
    # Activate the command and return the result.
    return `$commandString`;
}


=head3 CommaFormat

    my $formatted = Tracer::CommaFormat($number);

Insert commas into a number.

=over 4

=item number

A sequence of digits.

=item RETURN

Returns the same digits with commas strategically inserted.

=back

=cut

sub CommaFormat {
    # Get the parameters.
    my ($number) = @_;
    # Pad the length up to a multiple of three.
    my $padded = "$number";
    $padded = " " . $padded while length($padded) % 3 != 0;
    # This is a fancy PERL trick. The parentheses in the SPLIT pattern
    # cause the delimiters to be included in the output stream. The
    # GREP removes the empty strings in between the delimiters.
    my $retVal = join(",", grep { $_ ne '' } split(/(...)/, $padded));
    # Clean out the spaces.
    $retVal =~ s/ //g;
    # Return the result.
    return $retVal;
}


=head3 CompareLists

    my ($inserted, $deleted) = Tracer::CompareLists(\@newList, \@oldList, $keyIndex);

Compare two lists of tuples, and return a hash analyzing the differences. The lists
are presumed to be sorted alphabetically by the value in the $keyIndex column.
The return value contains a list of items that are only in the new list
(inserted) and only in the old list (deleted).

=over 4

=item newList

Reference to a list of new tuples.

=item oldList

Reference to a list of old tuples.

=item keyIndex (optional)

Index into each tuple of its key field. The default is 0.

=item RETURN

Returns a 2-tuple consisting of a reference to the list of items that are only in the new
list (inserted) followed by a reference to the list of items that are only in the old
list (deleted).

=back

=cut

sub CompareLists {
    # Get the parameters.
    my ($newList, $oldList, $keyIndex) = @_;
    if (! defined $keyIndex) {
        $keyIndex = 0;
    }
    # Declare the return variables.
    my ($inserted, $deleted) = ([], []);
    # Loop through the two lists simultaneously.
    my ($newI, $oldI) = (0, 0);
    my ($newN, $oldN) = (scalar @{$newList}, scalar @{$oldList});
    while ($newI < $newN || $oldI < $oldN) {
        # Get the current object in each list. Note that if one
        # of the lists is past the end, we'll get undef.
        my $newItem = $newList->[$newI];
        my $oldItem = $oldList->[$oldI];
        if (! defined($newItem) || defined($oldItem) && $newItem->[$keyIndex] gt $oldItem->[$keyIndex]) {
            # The old item is not in the new list, so mark it deleted.
            push @{$deleted}, $oldItem;
            $oldI++;
        } elsif (! defined($oldItem) || $oldItem->[$keyIndex] gt $newItem->[$keyIndex]) {
            # The new item is not in the old list, so mark it inserted.
            push @{$inserted}, $newItem;
            $newI++;
        } else {
            # The item is in both lists, so push forward.
            $oldI++;
            $newI++;
        }
    }
    # Return the result.
    return ($inserted, $deleted);
}

=head3 Cmp

    my $cmp = Tracer::Cmp($a, $b);

This method performs a universal sort comparison. Each value coming in is
separated into a text parts and number parts. The text
part is string compared, and if both parts are equal, then the number
parts are compared numerically. A stream of just numbers or a stream of
just strings will sort correctly, and a mixed stream will sort with the
numbers first. Strings with a label and a number will sort in the
expected manner instead of lexically. Undefined values sort last.

=over 4

=item a

First item to compare.

=item b

Second item to compare.

=item RETURN

Returns a negative number if the first item should sort first (is less), a positive
number if the first item should sort second (is greater), and a zero if the items are
equal.

=back

=cut

sub Cmp {
    # Get the parameters.
    my ($a, $b) = @_;
    # Declare the return value.
    my $retVal;
    # Check for nulls.
    if (! defined($a)) {
        $retVal = (! defined($b) ? 0 : -1);
    } elsif (! defined($b)) {
        $retVal = 1;
    } else {
        # Here we have two real values. Parse the two strings.
        my @aParsed = _Parse($a);
        my @bParsed = _Parse($b);
        # Loop through the first string.
        while (! $retVal && @aParsed) {
            # Extract the string parts.
            my $aPiece = shift(@aParsed);
            my $bPiece = shift(@bParsed) || '';
            # Extract the number parts.
            my $aNum = shift(@aParsed);
            my $bNum = shift(@bParsed) || 0;
            # Compare the string parts insensitively.
            $retVal = (lc($aPiece) cmp lc($bPiece));
            # If they're equal, compare them sensitively.
            if (! $retVal) {
                $retVal = ($aPiece cmp $bPiece);
                # If they're STILL equal, compare the number parts.
                if (! $retVal) {
                    $retVal = $aNum <=> $bNum;
                }
            }
        }
    }
    # Return the result.
    return $retVal;
}

# This method parses an input string into a string parts alternating with
# number parts.
sub _Parse {
    # Get the incoming string.
    my ($string) = @_;
    # The pieces will be put in here.
    my @retVal;
    # Loop through as many alpha/num sets as we can.
    while ($string =~ /^(\D*)(\d+)(.*)/) {
        # Push the alpha and number parts into the return string.
        push @retVal, $1, $2;
        # Save the residual.
        $string = $3;
    }
    # If there's still stuff left, add it to the end with a trailing
    # zero.
    if ($string) {
        push @retVal, $string, 0;
    }
    # Return the list.
    return @retVal;
}

=head3 ListEQ

    my $flag = Tracer::ListEQ(\@a, \@b);

Return TRUE if the specified lists contain the same strings in the same
order, else FALSE.

=over 4

=item a

Reference to the first list.

=item b

Reference to the second list.

=item RETURN

Returns TRUE if the two parameters are identical string lists, else FALSE.

=back

=cut

sub ListEQ {
    # Get the parameters.
    my ($a, $b) = @_;
    # Declare the return variable. Start by checking the lengths.
    my $n = scalar(@$a);
    my $retVal = ($n == scalar(@$b));
    # Now compare the list elements.
    for (my $i = 0; $retVal && $i < $n; $i++) {
        $retVal = ($a->[$i] eq $b->[$i]);
    }
    # Return the result.
    return $retVal;
}

=head3 Clean

    my $cleaned = Tracer::Clean($string);

Clean up a string for HTML display. This not only converts special
characters to HTML entity names, it also removes control characters.

=over 4

=item string

String to convert.

=item RETURN

Returns the input string with anything that might disrupt an HTML literal removed. An
undefined value will be converted to an empty string.

=back

=cut

sub Clean {
    # Get the parameters.
    my ($string) = @_;
    # Declare the return variable.
    my $retVal = "";
    # Only proceed if the value exists.
    if (defined $string) {
        # Get the string.
        $retVal = $string;
        # Clean the control characters.
        $retVal =~ tr/\x00-\x1F/?/;
        # Escape the rest.
        $retVal = CGI::escapeHTML($retVal);
    }
    # Return the result.
    return $retVal;
}

=head3 SortByValue

    my @keys = Tracer::SortByValue(\%hash);

Get a list of hash table keys sorted by hash table values.

=over 4

=item hash

Hash reference whose keys are to be extracted.

=item RETURN

Returns a list of the hash keys, ordered so that the corresponding hash values
are in alphabetical sequence.

=back

=cut

sub SortByValue {
    # Get the parameters.
    my ($hash) = @_;
    # Sort the hash's keys using the values.
    my @retVal = sort { Cmp($hash->{$a}, $hash->{$b}) } keys %$hash;
    # Return the result.
    return @retVal;
}

=head3 GetSet

    my $value = Tracer::GetSet($object, $name => $newValue);

Get or set the value of an object field. The object is treated as an
ordinary hash reference. If a new value is specified, it is stored in the
hash under the specified name and then returned. If no new value is
specified, the current value is returned.

=over 4

=item object

Reference to the hash that is to be interrogated or updated.

=item name

Name of the field. This is the hash key.

=item newValue (optional)

New value to be stored in the field. If no new value is specified, the current
value of the field is returned.

=item RETURN

Returns the value of the named field in the specified hash.

=back

=cut

sub GetSet {
    # Get the parameters.
    my ($object, $name, $newValue) = @_;
    # Is a new value specified?
    if (defined $newValue) {
        # Yes, so store it.
        $object->{$name} = $newValue;
    }
    # Return the result.
    return $object->{$name};
}

1;
