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
    use Getopt::Long::Descriptive;

=head1 Apache VHosts Fixup Script

    ApacheFix hostFileName

This script adds support for the C<fig.localhost> site to the Apache configuration.
It is used in Eclipse configurations to provide quick access to documentation tools.

=head2 Parameters

The single positional parameter is the name of the VHosts file. On a Mac, this is
C</etc/apache2/extra/httpd=vhosts.conf>. On Windows, this will depend on which
Apache installation you are using. For example, on a standard XAMPP installation,
it would be C</xampp/apache/conf/extra/httpd-vhosts.conf>.

The command-line options are as follows:

=over 4

=item mac

If this is a mac, the location of the main Apache configuration file. This is usually
C</private/etc/apache2/httpd.conf>.

=back

On a Mac, you must have root privileges to run this script (which is why it is no
longer a part of L<Config.pl>).

=cut

    $| = 1; # Prevent buffering on STDOUT.
    # Get the command-line options.
    my ($opt, $usage) = describe_options('%o %c vhostsFile',
    		['mac=s', 'main Apache configuration file (if this is a Mac)']);
    my $fileName = $ARGV[0];
    if (! $fileName) {
        die "The VHOSTS file name is required.";
    }
    # Determine the operating system.
    my $winMode = ($^O =~ /Win/ ? 1 : 0);
    # We'll put the file lines in here, omitting any existing SEEDtk section.
    my @lines;
    # Set up the VHOSTS file. Is there a previous copy?
    if (-f $fileName) {
	    # Open the configuration file for input.
	    open(my $ih, "<$fileName") || die "Could not open configuration file $fileName: $!";
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
    }
    # Now the file lines we want to keep from the old file (if any) are
    # in @lines. Open the file for output.
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
    # Check for a Mac requirement.
    if ($opt->mac) {
    	# Here we must update the main config file to enable virtual hosting.
		my $confFile = $opt->mac;
    	open(my $ih, "<$confFile") || die "Could not open $confFile: $!";
		# We'll accumulate output lines in here.
		@lines = ();
		# This will be set to the number of lines updated.
		my $count = 0;
    	while (! eof $ih) {
    		my $line = <$ih>;
    		if ($line =~ /^(\s*)#(\w+\s.+vhost.*)/ ||
    			$line =~ /^(\s*)#(LoadModule\s+cgi_module.+)/ ||
    			$line =~ /^(\s*)#(AddHandler\s+cgi.+)/) {
    			# Here we want to uncomment the line.
    			push @lines, "$1$2\n";
    			$count++;
    		} else {
    			# Here the line should be kept unchanged.
    			push @lines, $line;
    		}
    	}
    	# Close the file.
    	close $ih;
    	# Do we need to update it?
    	if ($count) {
			# Yes. Open it for output.
	    	undef $oh;
	    	open($oh, ">$confFile") || die "Could not open $confFile for output: $!";
			# Write out the updated lines.
	    	for my $line (@lines) {
	    		print $oh $line;
	    	}
	    	close $oh;
	    	print "$count lines updated in $confFile.\n";
    	}
    }

