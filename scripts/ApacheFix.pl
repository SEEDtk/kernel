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
    use FIG_Config;
    use Env;
    use File::Spec;

=head1 Apache VHosts Fixup Script

    ApacheFix hostFileName

This script adds support for the C<fig.localhost> site to the Apache configuration.
It is used in Eclipse configurations to provide quick access to documentation tools.

=head2 Parameters

The single positional parameter is the name of the VHosts file. On a Mac, this is
C</etc/apache2/extra/httpd=vhosts.conf>. On Windows, this will depend on which
Apache installation you are using. For example, on a standard XAMPP installation,
it would be C</xampp/apache/conf/extra/httpd-vhosts.conf>.

On a Mac, you must have root privileges to run this script (which is why it is no
longer a part of L<Config.pl>).

=cut

    $| = 1; # Prevent buffering on STDOUT.
    my $fileName = $ARGV[0];
    if (! $fileName) {
        die "The VHOSTS file name is required.";
    } elsif (! -f $fileName) {
        die "$fileName does not exist.";
    }
    # Determine the operating system.
    my $winMode = ($^O =~ /Win/ ? 1 : 0);
    # Set up the VHOSTS file.
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
    # Now we add our new stuff. First, get the name of the web directory.
    my $webdir = File::Spec->rel2abs($FIG_Config::web_dir);
    # Rel2Abs added a drive letter if we needed it, but we must fix the Windows
    # backslash craziness. Apache requires forward slashes.
    $webdir =~ tr/\\/\//;
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
    if ($winMode) {
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
    print "VHOSTS file updated.\n";

