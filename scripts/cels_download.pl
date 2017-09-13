#!/usr/bin/env perl
#
# Copyright (c) 2003-2015 University of Chicago and Fellowship
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
use warnings;
use File::Copy::Recursive;

my $dir = "~parrello/SEEDtk/Data/Cancer/CellLines";
File::Copy::Recursive::pathmk($dir) || die "Could not create directory: $!";
chdir($dir);
my @curlOpts = ('--remote-name', '--silent', '--show-error');
my $url = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3610/";
for (my $i = 1; $i <= 25; $i++) {
    my $file = "E-MTAB-3610.raw.$i.zip";
    print "Downloading $file.\n";
    my $rc = system('curl', @curlOpts, "$url$file");
    die "Error code $rc downloading file $i." if $rc;
    print "Unpacking $file.\n";
    $rc = system('unzip', "$file");
    die "Error code $rc unpacking $file." if $rc;
    print "Deleting $file.\n";
    kill $file;
}
