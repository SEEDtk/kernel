#!/usr/bin/env perl

#
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
#


    use strict;
    use File::Basename;
    use File::Spec;
    use File::Copy;
    use File::Path;
    use Getopt::Long::Descriptive;

    # We need to look inside the FIG_Config even though it is loaded at
    # run-time, so we will get lots of warnings about one-time variables.
    no warnings qw(once);

    ## THESE TWO CONSTANTS DEFINE THE PUBLIC SCRIPT and PM LIBRARIES.
    use constant SCRIPTS => qw(utils);
    use constant LIBS => qw(config ERDB kernel);

=head1 Generate SEEDtk Configuration Files

    Config.pl [ options ] dataDirectory webDirectory

This method generates (or re-generates) the L<FIG_Config> and B<UConfig.sh> files for a
SEEDtk environment.

=head2 Parameters

The positional parameters are the location of the data folder and the location
of the web folder (see L<ReadMe> for more information about SEEDtk folders). If a
B<FIG_Config> file already exists, this information is not needed-- the existing values
will be used.

The command-line options are as follows.

=over 4

=item fc

If specified, the name of the B<FIG_Config> file for the output. If the name is specified
without a path, it will be put in the main C<config> folder. If C<off>, no B<FIG_Config>
file will be written.

=item uc

If specified, the name of the UConfig file for the output. If the name is specified
without a path, it will be put in the parent of the source directory.  If C<off>, no UConfig file
will be written. If C<sys>, the changes will be made directly to the environment (via
the registry). This last is only possible under Windows.

=item winmode

If C<1>, the system will be configured for Windows; if C<0>, the system will be configured for Unix.
If unspecified, the current operating system will be interrogated.

=item clear

If specified, the current B<FIG_Config> values will be ignored, and the configuration information will
be generated from scratch.

=item links

If specified, a prototype C<Links.html> file will be generated in the web directory if one does not
already exist.

=item apache

If specified, the name of an Apache web configuration file (usually the B<extras/httpd-vhosts.conf>),
which will be updated to run the SEEDtk development testing server.

=item dirs

If specified, the default data and web subdirectories will be set up.

=item home

If specified, a link named C<SEEDtk> will be placed in this directory that points to the web root
directory. The specified directory whould be the user's home web directory, and this allows the
SEEDtk web to be accessed from it. Only Unix systems should use this feature.

=item pfix

If C<1>, the system requires a PERL path fixup. The contents of C<PerlPath.sh> will be added
to the B<UConfig> file. This defaults to C<1> for Unix and C<0> for Windows.

=item dna

If specified, the location for the DNA repository. If you are using a shared database, then
this insures you are using the same repository as everyone else.

=back

=head2 Notes for Programmers

To add a new L<FIG_Config> parameter, simply add a call to L<Env/WriteParam> to the
L</WriteAllParams> method.

To change the B<UConfig> file, simply add a call to L<Env/WriteConfig> to the
L</WriteAllConfigs> method.


=cut

    $| = 1; # Prevent buffering on STDOUT.
    print "Retrieving current configuration.\n";
    # Get the directory this script is running it.
    my $base_dir;
    if ($ENV{KB_TOP})
    {
        $base_dir = $ENV{KB_TOP};
    } else {
	$base_dir = dirname(File::Spec->rel2abs(__FILE__));

        # Fix Windows slash craziness.
        $base_dir =~ tr/\\/\//;
        # Chop off the folder name to get the source root.
        $base_dir =~ s/\/\w+$//;
        # Load in the "Env" module.
        unshift @INC, "$base_dir/kernel";
    }
    
    require Env;
    # Determine the operating system.
    my $winMode = ($^O =~ /Win/ ? 1 : 0);
    # Analyze the command line.
    my ($opt, $usage) = describe_options('%o %c dataRootDirectory webRootDirectory',
            ["winmode|w=i", "\"1\" for Windows, \"0\" for Unix", { default => $winMode }],
            ["clear|c", "ignore current configuration values"],
            ["fc=s", "name of a file to use for the FIG_Config output, or \"off\" to turn off FIG_Config output",
                    { default => "$base_dir/config/FIG_Config.pm" }],
            ["uc=s", "name of a file to use for the UConfig output, \"off\" to turn off UConfig output, or \"sys\" to write directly to the environment",
                    { default => "UConfig.sh" }],
            ["apache=s", "location of the Apache configuration files for the testing server"],
            ["dirs", "verify default subdirectories exist"],
            ["pfix=i", "perform PERL path fixup", { default => (1 - $winMode) }],
            ["dna=s", "location of the DNA repository (if other than local)"],
            ["home=s", "location of the home web directory for the current user; if specified, a link will be placed to the web directory"],
            ["links", "generate Links.html file"],
            );
    # Validate the web directory options.
    my $apache = $opt->apache;
    my $homeDir = $opt->home;
    if ($apache) {
        if (! -f $apache) {
            die "Apache config file $apache not found."
        } elsif ($homeDir) {
            die "Cannot have both the --apache and --home options."
        }
    } elsif ($homeDir) {
        if (! -d $homeDir) {
            die "Home directory $homeDir not found.";
        } elsif ($winMode) {
            die "Home directory option not valid for Windows.\n";
        }
    }
    # Check for mutually exclusive options.
    if ($opt->dirs && $winMode != $opt->winmode) {
        die "Option --dirs prohibited when targeting cross-platform.";
    } elsif ($homeDir && $opt->winmode) {
        die "Option --home prohibited when targeting Windows.";
    } elsif ($opt->uc eq 'sys' && $opt->pfix) {
        die "Cannot use PERL fixup (pfix) when doing a registry update.";
    }
    print "Analyzing directories.\n";
    # The root directories will be put in here.
    my ($dataRootDir, $webRootDir) = ('', '');
    # Get the name of the real FIG_Config file (not the output file,
    # if one was specified, the real one).
    my $fig_config_name = "$base_dir/config/FIG_Config.pm";
    # Now we want to get the current environment. If the CLEAR option is
    # specified or there is no file present, we stay blank; otherwise, we
    # load the existing FIG_Config.
    if (! $opt->clear && -f $fig_config_name) {
        RunFigConfig($fig_config_name);
    }
    # Make sure we have the data directory if there is no data root
    # in the command-line parameters.
    if (! defined $FIG_Config::shrub_dir) {
        $dataRootDir = FixPath($ARGV[0]);
        if (! defined $dataRootDir) {
            die "A data root directory is required if no current value exists in FIG_Config.";
        } elsif (! -d $dataRootDir) {
            die "The specified data root directory $dataRootDir was not found.";
        }
        # Are we setting up default data directories?
        if ($opt->dirs) {
            # Yes. Insure we have the data paths.
            BuildPaths($opt->winmode, Data => $dataRootDir, qw(Inputs Inputs/GenomeData Inputs/SubSystemData LoadFiles));
            # Are we using a local DNA repository?
            if (! $opt->dna) {
                # Yes. Build that, too.
                BuildPaths($opt->winmode, Data => $dataRootDir, qw(DnaRepo));
            }
        }
    }
    # Make sure we have the web directory if there is no web root in
    # the command-line parameters.
    if (! defined $FIG_Config::web_dir) {
        $webRootDir = FixPath($ARGV[1]);
        if (! defined $webRootDir) {
            die "A web root directory is required if no current value exists in FIG_Config.";
        } elsif (! -d $webRootDir) {
            die "The specified web root directory $webRootDir was not found.";
        }
        # Are we setting up default web directories?
        if ($opt->dirs) {
            # Yes. Insure we have the web paths.
            BuildPaths($opt->winmode, Web => $webRootDir, qw(img Tmp logs));
        }
    }
    #If the FIG_Config write has NOT been turned off, then write the FIG_Config.
    if ($opt->fc eq 'off') {
        print "FIG_Config output suppressed.\n";
    } else {
        # Compute the FIG_Config file name.
        my $outputName = $opt->fc;
        # Fix the slash craziness for Windows.
        $outputName =~ tr/\\/\//;
        # If the name is pathless, put it in the config directory.
        if ($outputName !~ /\//) {
            $outputName = "$base_dir/config/$outputName";
        }
        # If we are overwriting the real FIG_Config, back it up.
        if ($outputName eq $fig_config_name) {
            copy $fig_config_name, "$base_dir/config/FIG_Config_old.pm";
        }
        # Write the FIG_Config.
        print "Writing configuration to $outputName.\n";
        WriteAllParams($outputName, $base_dir, $dataRootDir, $webRootDir, $opt);
        # Execute it to get the latest variable values.
        print "Reading back new configuration.\n";
        RunFigConfig($outputName);
    }
    # Create the web configuration file. We need the key directories.
    my $sourcedir = $FIG_Config::source;
    my $webConfig = "$FIG_Config::web_dir/lib/Web_Config.pm";
    # Open the web configuration file for output.
    if (! open(my $oh, ">$webConfig")) {
        # Web system problems are considered warnings, not fatal errors.
        warn "Could not open web configuration file $webConfig: $!\n";
    } else {
        # Write the file.
        print $oh "\n";
        print $oh "    use lib\n";
        print $oh "        '" .	join("',\n        '", @FIG_Config::libs) . "';\n";
        print $oh "\n";
        print $oh "    use FIG_Config;\n";
        print $oh "\n";
        print $oh "1;\n";
        # Close the file.
        close $oh;
        print "Web configuration file $webConfig created.\n";
    }
    # If the UConfig write has NOT been turned off, then write the UConfig.
    if ($opt->uc eq 'off') {
        print "UConfig output suppressed.\n";
    } else {
        # Compute the output file name.
        my $ucFileName = $opt->uc;
        # If it is not a special name, we need to normalize it.
        if ($ucFileName ne 'sys') {
            # Fix the slash craziness for Windows.
            $ucFileName =~ tr/\\/\//;
            # If the name is pathless, put it in the source directory's parent.
            if ($ucFileName !~ /\//) {
                my $ucDir = $base_dir;
                $ucDir =~ s/\/\w+$//;
                $ucFileName = "$ucDir/$ucFileName";
            }
        }
        # Write the UConfig.
        WriteAllConfigs($ucFileName, $base_dir, $opt);
        print "UConfig file $ucFileName updated.\n";
        # Now we create the run_perl file in the web directory.  First, we
        # open the file for output.
        if (! open(my $oh, ">$FIG_Config::web_dir/run_perl.sh")) {
            # Web configuration problems are considered warnings, not fatal errors.
            warn "Could not update web run_perl file: $!";
        } else {
            # Are we using a PERL fixup?
            if ($opt->pfix) {
                # Yes. Copy in the PERL path fixup.
                CopyPerlFix($oh, $base_dir);
            }
            # Put in the command to execute PERL.
            print $oh 'exec perl "$@"' . "\n";
            # Close the output.
            close $oh;
            print "Execution helper run_perl.sh created.\n";
        }
        if (! $winMode && ! $opt->winmode) {
            # Here we are on a Unix system and not targeting Windows, so we need to
            # fix the execution permissions of all the script files.
            print "Fixing execution permissions.\n";
            # Fix the web directory permissions.
            FixPermissions($FIG_Config::web_dir, ".cgi", 0111);
            # Fix the permissions of the script directories.
            for my $scriptDir (@FIG_Config::scripts) {
                FixPermissions($scriptDir, ".pl", 0111);
            }
        }
    }
    # Check for an Apache Vhosts update request.
    if ($apache) {
        # Yes. Do the update.
        print "Updating Apache configuration file $apache.\n";
        SetupVHosts($apache, $opt->winmode);
    }
    # Check for the home directory symlink request.
    if ($homeDir) {
        my $linkDir = "$homeDir/SEEDtk";
        if (-d $linkDir) {
            print "Symbolic link directory $linkDir already exists.\n";
        } else {
            symlink($FIG_Config::web_dir, "$homeDir/SEEDtk");
            print "Symbolic link at $homeDir created for web hosting.\n";
        }
    }
    # Finally, check for the links file.
    if ($opt->links) {
        # Determine the output location for the links file.
        my $linksDest = "$FIG_Config::web_dir/Links.html";
        # Do we need to generate a links file?
        if (-f $linksDest) {
            # No need. We already have one.
            print "$linksDest file already exists-- not updated.\n";
        } else {
            # We don't have a links file yet.
            print "Generating new $linksDest.\n";
            # Find the source copy of the file.
            my $linksSrc = "$base_dir/utils/Links.html";
            # Copy it to the destination.
            copy $linksSrc, $linksDest;
            print "$linksDest file created.\n";
        }
    }
    print "All done.\n";

=head2 Internal Subroutines

=head3 RunFigConfig

    RunFigConfig($fileName);

Execute the L<FIG_Config> module. This uses the PERL C<do> function, which
unlike C<require> can execute a module more than once, but requires error
checking. The error checking is done by this method.

=over 4

=item fileName

The name of the B<FIG_Config> file to load.

=back

=cut

sub RunFigConfig {
    # Get the parameters.
    my ($fileName) = @_;
    # Execute the FIG_Config;
    do $fileName;
    if ($@) {
        # An error occurred compiling the module.
        die "Error compiling FIG_Config: $@";
    } elsif ($!) {
        # An error occurred reading the module.
        die "Error reading FIG_Config: $!";
    }
}

=head3 WriteAllParams

    WriteAllParams($fig_config_name, $base_dir, $dataRootDir, $opt);

Write out the B<FIG_Config> file to the specified location. This method
is mostly calls to the L</WriteParam> method, which provides a concise
way of writing parameters to the file and checking for pre-existing
values. It is presumed that L</RunFigConfig> has been executed first so
that the existing values are known.

=over 4

=item fig_config_name

File name for the B<FIG_Config> file. The parameter code will be written to
this file.

=item base_dir

Location of the base directory for the source code.

=item dataRootDir

Location of the base directory for the data and web files.

=item opt

Command-line options object.

=back

=cut

sub WriteAllParams {
    # Get the parameters.
    my ($fig_config_name, $base_dir, $dataRootDir, $webRootDir, $opt) = @_;
    # Open the FIG_Config for output.
    open(my $oh, ">$fig_config_name") || die "Could not open $fig_config_name: $!";
    # Write the initial lines.
    print $oh "package FIG_Config;\n";
    Env::WriteLines($oh,
        "",
        "## WHEN YOU ADD ITEMS TO THIS FILE, BE SURE TO UPDATE source/utils/Config.pl.",
        "## All paths should be absolute, not relative.",
        "");
    # Write each parameter.
    Env::WriteParam($oh, 'root directory of the local web server', web_dir => $webRootDir);
    Env::WriteParam($oh, 'directory for temporary files', temp => "$webRootDir/Tmp");
    Env::WriteParam($oh, 'URL for the directory of temporary files', temp_url => 'http://localhost/Tmp');
    Env::WriteParam($oh, 'TRUE for windows mode', win_mode => ($opt->winmode ? 1 : 0));
    Env::WriteParam($oh, 'source code root directory', source => $base_dir);
    Env::WriteParam($oh, 'BASH prompt string', prompt => '\h-\W');
    ## Put new non-Shrub parameters here.
    # These next parameters are lists, so we have to build them manually.
    Env::WriteLines($oh, "", "# base names of the private script folders",
            'our @pscripts = qw(' . join(" ", @FIG_Config::pscripts) . ');',
            "", "# base names of the private PERL library folders",
            'our @plibs = qw(' . join(" ", @FIG_Config::plibs) . ');');
    # These next parameters use code, so are not subject to the usual Env::WriteParam logic.
    Env::WriteLines($oh, "", "# list of script directories",
            'our @scripts = map { "$source/$_" } (qw(' . join(" ", SCRIPTS) . '), @pscripts);',
            "",  "# list of PERL libraries",
            'our @libs = map { "$source/$_" } (qw(' . join(" ", LIBS) . '), @plibs);');
    # Now comes the Shrub configuration section.
    Env::WriteLines($oh, "", "", "# SHRUB CONFIGURATION", "");
    Env::WriteParam($oh, 'root directory for Shrub data files (should have subdirectories "Inputs" (optional) and "LoadFiles" (required))',
            shrub_dir => "$dataRootDir");
    Env::WriteParam($oh, 'full name of the Shrub DBD XML file', shrub_dbd => "$base_dir/ERDB/ShrubDBD.xml");
    Env::WriteParam($oh, 'Shrub database signon info (name/password)', userData => "seed/");
    Env::WriteParam($oh, 'name of the Shrub database (empty string to use the default)', shrubDB => "");
    Env::WriteParam($oh, 'TRUE if we should create indexes before a table load (generally TRUE for MySQL, FALSE for PostGres)',
            preIndex => 1);
    Env::WriteParam($oh, 'default DBMS (currently only "mysql" works for sure)', dbms => "mysql");
    Env::WriteParam($oh, 'database access port', dbport => 3306);
    Env::WriteParam($oh, 'TRUE if we are using an old version of MySQL (legacy parameter; may go away)', mysql_v3 => 0);
    Env::WriteParam($oh, 'default MySQL storage engine', default_mysql_engine => "InnoDB");
    Env::WriteParam($oh, 'database host server (empty string to use the default)', dbhost => "");
    Env::WriteParam($oh, 'TRUE to turn off size estimates during table creation-- should be FALSE for MyISAM',
            disable_dbkernel_size_estimates => 1);
    Env::WriteParam($oh, 'mode for LOAD TABLE INFILE statements, empty string is OK except in special cases (legacy parameter; may go away)',
            load_mode => '');
    # Now comes the DNA repository. There are two cases-- a local repository or a global one. Check the type.
    my $dnaRepo = $opt->dna;
    if (! $dnaRepo) {
        # Here we have the local repository.
        $dnaRepo = "$dataRootDir/DnaRepo";
    }
    Env::WriteParam($oh, 'location of the DNA repository', shrub_dna => $dnaRepo);
    ## Put new Shrub parameters here.
    # Write the trailer.
    print $oh "\n1;\n";
    # Close the output file.
    close $oh;
}

=head3 WriteAllConfigs

    WriteAllConfigs($fileName, $base_dir, $opt);

Write out the B<UConfig.sh> file. This file is used in Unix systems to
set up environment variables for PERL includes and execution paths. It
presumes the B<FIG_Config> file has been updated and L</RunFigConfig> has
been called to load its variables.

=over 4

=item fileName

Name of the output file. If C<sys>, then the updates will be made directly to
the system environment. This is only possible if the script is being run in
administrator mode in a Windows system.

=item base_dir

Name of the base directory for the source.

=item opt

Command-line options object.

=back

=cut

    use constant REGKEY => 'SYSTEM/CurrentControlSet/Control/Session Manager/Environment';

sub WriteAllConfigs {
    # Get the parameters.
    my ($fileName, $base_dir, $opt) = @_;
    # This will be the delimiter used for path merging. This is
    # a colon for Unix, semi-colon for Windows.
    my $delim = ($opt->winmode ? ';' : ':');
    # This will be the output file handle. We'll set it to the
    # environment cluster if we are writing to the registry.
    my $oh;
    if ($opt->uc eq 'sys') {
        # Here we are updating the registry.
        $oh = Env->GetRegKey(HKLM => REGKEY);
        print "Writing environment changes to registry.\n";
    } else {
        # Here we are writing to a file.
        open($oh, ">$fileName") || die "Could not open UConfig file $fileName: $!";
        print "Writing environment changes to $fileName.\n";
        # Check for a PERL fix requirement.
        if ($opt->pfix) {
            CopyPerlFix($oh, $base_dir);
        }
    }
    # Compute the script paths.
    my $paths = Env::BuildPathList($opt->winmode, $delim, @FIG_Config::scripts);
    Env::WriteConfig($oh, "Add SEEDtk scripts to the execution path.", PATH => $paths, merge => $delim, expanded => 1);
    # Set the PERL libraries.
    my $libs = Env::BuildPathList($opt->winmode, $delim, @FIG_Config::libs);
    Env::WriteConfig($oh, "Add SEEDtk libraries to the PERL library path.", PERL5LIB => $libs, merge => $delim);
    ## Put new configuration parameters here.
    # This next section is for shell scripts only.
    if (ref $oh eq 'GLOB') {
        # Change to the SEEDtk directory.
        my $rootDir = $base_dir;
        $rootDir =~ s/\/source$//;
        print $oh "cd $rootDir\n";
        # Fix up the prompt
        print $oh "PS1=\"$FIG_Config::prompt\$ \"\n";
    }
    # The file (or registry key) in $oh will close automatically when we go out of scope.
}

=head3 CopyPerlFix

    CopyPerlFix($oh, $base_dir);

Copy the B<PerlPath.sh> file into the specified output file. This puts a PERL directory correction
into the path so that the correct PERL version is run.

=over 4

=item oh

Output handle to which the PERL fixup should be written.

=item base_dir

Root directory for the source tree. B<PerlPath.sh> will be in its B<config> subdirectory.

=back

=cut

sub CopyPerlFix {
    # Get the parameters.
    my ($oh, $base_dir) = @_;
     # Start with a comment.
       print $oh "# Fix PERL execution path.\n";
       # Get the PERL fix file.
       open(my $ih, "<$base_dir/config/PerlPath.sh") || die "Could not open PerlPath.sh: $!";
    # Spool it to the environment script.
    print $oh (<$ih>);
    # Add a spacer.
    print $oh "\n";
}


=head3 SetupVHosts

    SetupVHosts($fileName, $winmode);

This method updates an Apache vhosts file to enable the fig.localhost virtual local
server for SEEDtk development. It replaces or adds a specially-marked section that
defines the server root and its characteristics.

=over 4

=item fileName

Name of the B<vhosts.conf> file to udpate.

=item winmode

TRUE if the server is on Windows, else FALSE.

=back

=cut

sub SetupVHosts {
    # Get the parameters.
    my ($fileName, $winmode) = @_;
    # Open the configuration file for input.
    open(my $ih, "<$fileName") || die "Could not open configuration file $fileName: $!";
    # We'll put the file lines in here, omitting any existing SEEDtk section.
    my @lines;
    my $skipping;
    while (! eof $ih) {
        my $line = <$ih>;
        # Are we in the SEEDtk section?
        if ($skipping) {
            # Yes. Check for an end marker.
            if ($line =~ /^## END SEEDtk SECTION/) {
                # Found it. Stop skipping.
                $skipping = 0;
            }
        } else {
            # No. Check for a begin marker.
            if ($line =~ /^## BEGIN SEEDtk SECTION/) {
                # Found it. Start skipping.
                $skipping = 1;
            } else {
                # Not a marker. Save the line.
                push @lines, $line;
            }
        }
    }
    # Close the file.
    close $ih;
    # Open it again for output.
    open(my $oh, ">$fileName") || die "Could not open configuration file $fileName: $!";
    # Unspool the lines from the old file.
    for my $line (@lines) {
        print $oh $line;
    }
    # Now we add our new stuff. First, get the name of the web and source directories.
    # The BuildPathList normalizes the paths according to the target environment.
    my $paths = Env::BuildPathList($winmode, ";", $FIG_Config::web_dir, $FIG_Config::source);
    # Fix the Windows backslash craziness. Apache requires forward slashes.
    $paths =~ tr/\\/\//;
    # Extract the individual directories.
    my ($webdir, $sourcedir) = split /;/, $paths;
    # Write the start marker.
    print $oh "## BEGIN SEEDtk SECTION\n";
    # Declare the root directory for the virtual host.
    print $oh "<Directory \"$webdir\">\n";
    print $oh "    Options Indexes FollowSymLinks ExecCGI\n";
    print $oh "    AllowOverride None\n";
    print $oh "    Require all granted\n";
    print $oh "</Directory>\n";
    print $oh "\n";
    # Configure the virtual host itself.
    print $oh "<VirtualHost *:80>\n";
    # Declare the URL and file location of the root directory.
    print $oh "    DocumentRoot \"$webdir\"\n";
    print $oh "    ServerName fig.localhost\n";
    # If this is Windows, set up the registry for CGI execution.
    if ($winmode) {
        print $oh "    ScriptInterpreterSource Registry\n";
    }
    # Define the local logs.
    print $oh "    ErrorLog \"$webdir/logs/error.log\"\n";
    print $oh "    CustomLog \"$webdir/logs/access.log\" common\n";
    # Set up the default files for each directory to the usual suspects.
    print $oh "    DirectoryIndex index.cgi index.html index.htm\n";
    # Finish the host definition.
    print $oh "</VirtualHost>\n";
    # Write the end marker.
    print $oh "## END SEEDtk SECTION\n";
    # Close the output file.
    close $oh;
}

=head3 BuildPaths

    BuildPaths($winmode, $label => $rootDir, @subdirs);

Create the desired subdirectories for the specified root directory. The
type of root directory is provided as a label for status messages. On
Unix systems, a C<chmod> will be performed to fix the privileges.

=over 4

=item winmode

TRUE for a Windows system, FALSE for a Unix system.

=item label

Label describing the type of directory being created.

=item rootDir

Root directory path. All new paths created will be under this one.

=item subdirs

List of path names, relative to the specified root directory, that we
must insure exist.

=back

=cut

sub BuildPaths {
    # Get the parameters.
    my ($winmode, $label, $rootDir, @subdirs) = @_;
    # Loop through the new paths.
    for my $path (@subdirs) {
        my $newPath = "$rootDir/$path";
        # Check to see if the directory is already there.
        if (! -d $newPath) {
            # No, we must create it.
            File::Path::make_path($newPath);
            print "$label directory $newPath created.\n";
            # If this is Unix, fix the permissions.
            if (! $winmode) {
                chmod 0777, $newPath;
            }
        }
    }
}

=head3 FixPermissions

    FixPermissions($directory, $ext, $mask);

Add the specified permissions to all the files in a directory that match
a given file extension. The extension should include the separating
period. So, for example, C<.pl> would match all PERL script files.

=over 4

=item directory

Name of the directory whose files are to be changed.

=item ext

Suffix for the files to be changed. All file names that end in this string
will be updated.

=item mask

Mask to be merged into the file permission mask. All bits that match the 1 bits
in this mask will be turned on.

=back

=cut

sub FixPermissions {
    # Get the parameters.
    my ($directory, $ext, $mask) = @_;
    # Open the directory.
    opendir(my $dh, $directory) || die "Could not open $directory: $!";
    # Get all the matching files.
    my $len = length($ext);
    my @files = grep { substr($_, -$len, $len) eq $ext } readdir($dh);
    # Close the directory.
    closedir $dh;
    # Loop through the files.
    for my $file (@files) {
        # Compose the full file name.
        my $fileName = "$directory/$file";
        # Compute the new mode.
        my @finfo = stat $fileName;
        my $newMode = ($finfo[2] & 0777) | $mask;
        # Update the file.
        chmod $newMode, $fileName;
    }
}

=head3 FixPath

    my $absPath = FixPath($path);

Convert a path from relative to absolute, convert backslashes to slashes, and remove the drive letter.
This makes the path suitable for passing around in PERL.

=over 4

=item path

Path to convert.

=item RETURN

Returna an absolute, canonical version of the path.

=back

=cut

sub FixPath {
    # Get the parameters.
    my ($path) = @_;
    # Convert the path to a canonical absolute.
    my $retVal = File::Spec->rel2abs($path);
    # Convert backslashes to slashes.
    $retVal =~ tr/\\/\//;
    # Remove the drive letter (if any).
    $retVal =~ s/^\w://;
    # Return the result.
    return $retVal;
}
