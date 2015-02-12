package Env;

    use strict;
    use File::Spec;


=head1 Environment Modification Utilities

This package contains utility methods for modifying the environment. This includes
methods for writing assignment statements to files and changing the Windows
registry. The Windows stuff is only loaded if it's required. We
use L<Win32::Registry> rather than the newer B<TieRegistry> module because
the newer one does not always work. Note that you must be running as an
administrator if you are trying to update the registry.

The SEEDtk environment requires win-bash and the GnuWin32 tools in order
to run under Windows. This means the configuration file formats are
identical between the two environments.

=head2 Registry Methods

=head3 GetRegKey

    my $keyStructure = Environment->GetRegKey($root => $path);

Return a L<Win32::Registry> structure for manipulating the values of a
particular registry key.

=over 4

=item root

Root of the registry section of interest. Permissible values are C<HKLM>
(HKEY_LOCAL_MACHINE), C<HKCR> (HKEY_CLASSES_ROOT), and C<HKCU>
(HKEY_CURRENT_USER).

=item path

Path from the root to the registry key. Either slashes C</> or
backslashes C<\> can be used as delimiters between path segments.

=item RETURN

Returns a L<Win32::Registry> object for manipulating the key in question.

=back

=cut

sub GetRegKey {

    # Get the parameters.
    my ($class, $root, $path) = @_;
    # Insure we have access to the registry stuff.
    require Win32::Registry;
    # Get the root variable.
    my $rootKey;
    {   # We temporarily turn off warnings in case we don't
        # load the Win32 stuff, in which case the variables
        # below are not declared.
        no warnings 'once';
        if ($root eq 'HKLM') {
            $rootKey = $::HKEY_LOCAL_MACHINE;
        } elsif ($root eq 'HKCR') {
            $rootKey = $::HKEY_CLASSES_ROOT;
        } elsif ($root = 'HKCU') {
            $rootKey = $::HKEY_CURRENT_USER;
        }
    }
    # Normalize the path. We convert slashes to backslashes and insure there is one at the end.
    $path =~ tr/\//\\/;
    if (substr($path, -1) ne "\\") {
        $path .= "\\";
    }
    # Ask for the registry key.
    my $retVal;
    $rootKey->Open($path, $retVal);
    # Return it to the caller.
    return $retVal;
}

=head2 Environment Query Methods

=head3 GetFigConfigs

    my $varHash = Env::GetFigConfigs($fileName);

Return a hash describing the variables in the specified L<FIG_Config> file.
The file will be parsed for comments (one or more lines beginning with
C<#>) followed by assignment statements. The assignment statements will
be parsed to get the variable name and type (which can be scalar C<$> or
list C<@>). The return hash will contain the name, comments, and value of
each variable. List values will be converted to comma-separated lists.
(Based on the type of thing in FIG_Config, this is not a high-risk
approach.)

=over 4

=item fileName

Name of the FIG_Config file.

=item RETURN

Returns a reference to a hash, keyed by variable name. Each variable is
mapped to a reference to a two-element list consisting of the comment
followed by the value.

=back

=cut

sub GetFigConfigs {
    # Get the parameters.
    my ($fileName) = @_;
    # Declare the return variable.
    my %retVal;
    # Open the file for input.
    open(my $ih, "<$fileName") || die "Could not open $fileName: $!";
    # This list will accumulate the current comment lines.
    my @comments;
    # Loop through the file.
    while (! eof $ih) {
        my $line = <$ih>;
        # Determine the line type.
        if ($line =~ /^\s*#\s*(.*)/) {
            # Here we have a comment line. Accumulate it in the comment list.
            push @comments, $1;
        } elsif ($line =~ /^\s*our\s+([\$|\@])(\w+)\s*=/) {
            # Here we have an assignment. Save the variable name and type.
            my ($varType, $varName) = ($1, $2);
            # Create the comment.
            my $comment = join(" ", @comments);
            # Form the full name of the variable.
            my $fullName = $varType . 'FIG_Config::' . $varName;
            # Convert it into an expression for the variable value. We only
            # need to do this for lists.
            if ($varType eq '@') {
                $fullName = '"(" . join(", ", ' . $fullName . ') . ")"';
            }
            # Get the desired value.
            my $value = eval($fullName);
            if ($@) {
                die "Error evaluating $varName: $@";
            }
            # Store the information about this variable in the hash.
            $retVal{"$varType$varName"} = [$comment, $value];
        } else {
            # Here we have a separator line. We must clear the
            # comment list.
            @comments = ();
        }
    }
    # Return the result.
    return \%retVal;
}

=head3 GetScripts

    my $scriptHash = Env::GetScripts($dir);

Return a hash of the script files in each directory.

=over 4

=item dir

Name of the directory shows scripts are to be computed.

=item RETURN

Returns a reference to a hash that maps each script name (without the path) to the title from its
POD documentation.

=back

=cut

sub GetScripts {
    # Get the parameters.
    my ($dir) = @_;
    # Declear the return variable.
    my %retVal;
    # Open the directory.
    opendir(my $dh, $dir) || die "Could not open directory $dir: $!";
    # Find all the script files.
    my @files = grep { $_ =~ /^\w+\.pl$/ } readdir($dh);
    close $dh;
    # Loop through the file names.
    for my $file (@files) {
        # Open this file for input.
        if (! open(my $ih, "<$dir/$file")) {
            # Here the open failed. Store the error message.
            $retVal{$file} = "Error reading file: $!";
        } else {
            # We opened the file. Look for the heading comment.
            my $comment;
            while (! eof $ih && ! $comment) {
                my $line = <$ih>;
                if ($line =~ /^=head1\s+(.+)/) {
                    $comment = $1;
                }
            }
            # Put the comment (or undef if there was none) in the hash.
               $retVal{$file} = $comment;
        }
    }
    # Return the hash.
    return \%retVal;
}


=head2 Parameter Output Methods

=head3 BuildPathList

    my $pathString = Env::BuildPathList($winMode, $delim, @paths);

Build a string that contains a list of file paths. The file paths will be
normalized into the appropriate form based on the target operating system.

=over 4

=item winMode

TRUE if the target is Windows, else FALSE.

=item delim

Delimiter to use for separating the paths.

=item paths

List of paths to join into the output string.

=item RETURN

Returns a delimited list of the path strings.

=back

=cut

sub BuildPathList {
    # Get the parameters.
    my ($winMode, $delim, @paths) = @_;
    # The normalized paths will be put in here.
    my @normals;
    for my $path (@paths) {
        # Is this Windows?
        if ($winMode) {
            # Yes. Convert the slashes.
            $path =~ tr/\//\\/;
            # Insure we have a drive letter.
            if ($path !~ /^\w:/) {
                $path = "C:$path";
            }
        }
        # Save the modified path.
        push @normals, $path;
    }
    # Form the result string.
    my $retVal = join($delim, @normals);
    # Return it.
    return $retVal;
}

=head3 WriteLines

    Env::WriteLines($oh, @lines);

Write one or more lines of text to the specified output file. Nonblank
lines will be indented four spaces.

=over 4

=item oh

Open file handle for writing to the B<FIG_Config> file.

=item lines

List of text lines to write. They will all be indented four spaces and followed by new-lines.

=back

=cut

sub WriteLines {
    # Get the parameters.
    my ($oh, @lines) = @_;
    # Loop through the lines.
    for my $line (@lines) {
        # If the line is nonblank, pad it on the left.
        if ($line =~ /\S/) {
            $line = "    $line";
        }
        # Write the line.
        print $oh "$line\n";
    }
}


=head3 WriteParam

    Env::WriteParam($oh, $comment, $varName => $defaultValue);

Write a parameter value to the B<FIG_Config>. If the parameter already
has a value, that will be written. Otherwise, the default value will be
written. The value will be preceded by the comment.

=over 4

=item oh

Open output file handle for the B<FIG_Config> file.

=item comment

Comment to write out before the assignment statement.

=item varName

Name of the variable to assign.

=item defaultValue

Value to assign to the variable if it does not already have one.

=back

=cut

    # Table for escaping quoted strings.
    use constant BACKSLASH => { ( map { $_ => "\\$_" } ( '\\', '"', '$', '@' ) ),
                                ( map { chr($_) => "\\x" . unpack('H2', chr($_)) } (0..31, 128..255) ),
                                ( "\r" => "\\r", "\n" => "\\n", "\t" => "\\t")
                              };

sub WriteParam {
    # Get the parameters.
    my ($oh, $comment, $varName, $defaultValue) = @_;
    # We will put the desired variable value in here.
    my $value = eval("\$FIG_Config::$varName");
    if ($@) {
        die "Error checking value of FIG_Config::$varName";
    }
    if (! defined $value) {
        # There is no existing value, so use the default.
        $value = $defaultValue;
    }
    # Convert the value to PERL.
    if ($value !~ /^\d+$/) {
        # Here the value is not a number. Convert it to backslashed form and quote it.
        my @output = ('"');
        for my $ch (split //, $value) {
            push @output, (BACKSLASH->{$ch} // $ch);
        }
        $value = join("", @output, '"');
    }
    # Create the list of comment lines. First, we split the comment
    # into words.
    my @words = split " ", $comment;
    # Build the comment lines in here.
    my @comments = ("#");
    my $i = 0;
    for my $word (@words) {
        if (length($comments[$i]) + length($word) > 60) {
            $comments[++$i] = "# $word";
        } else {
            $comments[$i] .= " $word";
        }
    }

    # Write a blank line followed by the parameter's comments and value.
    WriteLines($oh, "", @comments, "our \$$varName = $value;");
}


=head3 WriteConfig

    Env::WriteConfig($oh, $comment, $varName, $value, %options);

Write a configuration parameter. This may involve output to a file or
an update to the registry environment. Because we use a bash emulator
in Windows, the file output syntax is very similar in both ebvironments.

=over 4

=item oh

Either an open file handle to which we should write the parameter or
a L<Win32::Registry> object for managing the Windows environment strings.

=item comment

A description of what the environment variable means.

=item varName

The name of the environment variable to set.

=item value

The value to give to the environment variable.

=item options

A hash of options that modify the write operation. The following keys
are supported.

=over 8

=item merge

If specified, the value is presumed to be a list that is merged into the
configuration. The value of the option is the list delimiter. In the
registry update case, this is an imperfect system, since if the
values to be merged change, old data could be left behind.

=item expanded

If TRUE, then the value will be stored as an expanded string if it is
stored in the registry. Expanded strings allow %-variable substitution.

=back

=back

=cut

sub WriteConfig {
    # Get the parameters.
    my ($oh, $comment, $varName, $value, %options) = @_;
    # Check for a merge delimiter.
    my $delim = $options{merge};
    # Determine whether this is a registry update or a write.
    if (ref $oh ne 'Win32::Registry') {
        # Here we have an output file.
        # Write the comment.
        print $oh "# $comment\n";
        # Is this a merge?
        if ($delim) {
            # Yes. Write a command to prefix the new information.
            print $oh "$varName=\"$value$delim\$$varName\"\n";
            print $oh "export $varName\n";
        } else {
            # No. Write a command to store the new information.
            print $oh "$varName=\"$value\"\n";
            print $oh "export $varName\n";
        }
    } else {
        # Here we have a registry update. There is no comment. We'll
        # put the target value in here.
        my $actualValue = $value;
        # If this is a merge, we need to parse the existing value and
        # merge the new values in;
        if ($delim) {
            # Parse the incoming value.
            my %items = map { $_ => 1 } split $delim, $value;
            # Get the current environment value and parse it as well,
            # removing values found in @items.
            my @current = grep { ! $items{$_} } split $delim, $ENV{$varName};
            # Form a new value.
            $actualValue = join($delim, $value, @current);
        }
        # Determine the type of string being stored. 2 = expanded string,
        # 1 = string.
        my $type = ($options{expanded} ? 2 : 1);
        # Store the new value in the registry. The 0 in the parameter list
        # is a weird reserved parameter that has no effect but is still
        #required.
        $oh->SetValueEx($varName, 0, $type, $actualValue);
    }
}


1;