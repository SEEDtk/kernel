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
use FIG_Config;
use ScriptUtils;
use Stats;
use LWP::UserAgent;
use XML::Simple;
use HTTP::Request;

=head1 Download Data from the European Nucleotide Archive

    ena_list.pl [ options ] outputFile

This script downloads data from the European Nucleotide Archive. It accepts as input a tab-delimited file containing
data about the samples in a project.

The output file will contain the sample ID, selected attributes, and the FASTQ file name for each sample.

=head2 Parameters

The single positional parameter is the name of the output file. The standard output will be used for tracing and status.

The command-line options are those found in L<ScriptUtils/ih_options>.

=cut

use constant ATTRIBUTES => { alzheimers => 1, autoimmune => 1, bmi_cat => 1, cardiovascular_disease => 1,
        depression_bipolar_schizophrenia => 1, diabetes => 1, epilepsy_or_seizure_disorder => 1,
        ibd => 1, ibs => 1, kidney_disease => 1, liver_disease => 1, lung_disease => 1, pku => 1 
        };
# Create the user agent.
my $ua = LWP::UserAgent->new();
# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outFile', ScriptUtils::ih_options(),
        );
# Open the input file.
my $ih = ScriptUtils::IH($opt->input);
# Open the output file.
my $oh;
my ($oFile) = @ARGV;
if (! $oFile) {
    die "No output file specified.";
} elsif (! open($oh, '>', $oFile)) {
    die "Could not open output file: $!";
}
# Create a header line.
my @hdrs = ('id', 'url', 'healthy', sort keys %{ATTRIBUTES()});
# Create the statistics object.
my $stats = Stats->new();
# Clear the header line.
my $line = <$ih>;
# Count the input lines.
my $count = 0;
# Track the responses.
my %response;
# Loop through the input file.
while (! eof $ih) {
    $line = <$ih>;
    $line =~ s/\r?\n$//;
    my @cols = split /\t/, $line;
    $stats->Add(lineIn => 1);
    # Get the sample ID and the FASTQ url.
    my $id = $cols[3];
    my $url = 'ftp://' . $cols[10];
    my $xurl = "http://www.ebi.ac.uk/ena/data/view/$id&display=xml&download=xml&filename=$id.xml";
    my $request = HTTP::Request->new(GET => $xurl);
    my $response = $ua->request($request);
    if ($response->code ne 200) {
        $stats->Add(errorResponse => 1);
        print "Error response for $id XML download: " . $response->message . "\n";
    } else {
        # Get the xml.
        my $xml = $response->content;
        my $xmlH = XMLin($xml);
        # Get the attribute list.
        my $attributesL = $xmlH->{SAMPLE}{SAMPLE_ATTRIBUTES}{SAMPLE_ATTRIBUTE};
        if (! $attributesL) {
            print "No attributes in $id.\n";
        } else {
            # This will be TRUE if the patient is diseased.
            my $healthy = "yes";
            # Save the useful attributes.
            my %attrH;
            for my $attribute (@$attributesL) {
                my $tag = $attribute->{TAG};
                my $value = $attribute->{VALUE};
                if (ATTRIBUTES->{$tag}) {
                    $attrH{$tag} = $value;
                    $stats->Add(attrFound => 1);
                    $response{$tag}{$value}++;
                    if ($value =~ /diagnosed/i) {
                        $healthy = "no";
                        $stats->Add(valUnhealthy => 1);
                    } elsif ($value =~ /^Un/) {
                        $stats->Add(valAmbiguous => 1);
                    } elsif ($value eq 'I do not have this condition') {
                        $stats->Add(valHealthy => 1);
                    } else {
                        $stats->Add(valOther => 1);
                    }
                } else {
                    $stats->Add(attrSkipped => 1);
                }
            }
            # Add the main data.
            $attrH{id} = $id;
            $attrH{url} = $url;
            $attrH{healthy} = $healthy;
            # Form the output line.
            my @out;
            for my $hdr (@hdrs) {
                push @out, $attrH{$hdr} // '';
            }
            print $oh join("\t", @out) . "\n";
            $stats->Add(lineOut => 1);
            $count++;
            if ($count % 100 == 0) {
                print "$count samples processed.\n";
            }
        }
    }
}
# Process the response counts.
for my $tag (sort keys %response) {
    my $values = $response{$tag};
    print "Responses for $tag.\n";
    for my $value (sort keys %$values) {
        my $count = $values->{$value} . " ";
        $count = " $count" while length($count) < 10;
        print "     $count $value\n";
    }
}
print "\n";
# All done. Print the stats.
print "All done.\n" . $stats->Show();

