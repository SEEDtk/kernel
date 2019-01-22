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
use WWW::Search;
use WWW::Search::PubMed;
use Template;

=head1 Retrieve Papers from PUBMED Based on Keywords

    pubmed_papers.pl [ options ] outFile keyword1 keyword2 ...

This script will perform a keyword search of the PUBMED database and produce an HTML page listed the ID, title, and abstract of each
paper returned.  An HTML page is produced so that the abstract can be made more readable.  The abstracts can then be examined for indications
of useful information and the ID link used to navigate to it.

=head2 Parameters

The first positional parameter is the name to give to the output file.  The remaining positional parameters are the keywords to use for the
search.  These will be passed to PUBMED conjunctively-- that is, strung together with the C<AND> operator.

The command-line options are as follows.

=over 4

=item template

The location of the web page template file.  The default is C<pubmed.tt> in the shared source library directory.

=cut

# Get the command-line parameters.
my $opt = ScriptUtils::Opts('outFile keyword1 keyword2 ... keywordN',
        ['template=s', 'web page template file', { default => "$FIG_Config::mod_base/RASTtk/lib/pubmed.tt" }],
        );
# Get the parameters.
my ($outFile, @keywords) = @ARGV;
# Verify the output file.
my $oh;
if (! $outFile) {
    die "No output file specified.";
} elsif (! open($oh, '>', $outFile)) {
    die "Could not open $outFile: $!";
}
# Verify the keywords.
if (! @keywords) {
    die "At least one keyword must be specified.";
}
# This list will contain the descriptors for each paper found.
my @papers;
# This hash will contain the variable values for the HTML template.
my %vars = (
    keyword_list => join(' ', @keywords),
    result_count => 0,
    papers => \@papers
);
# Perform the search.
print "Searching for papers with " . scalar(@keywords) . " keywords.\n";
my $s = WWW::Search->new('PubMed');
$s->native_query(join(' AND ', @keywords));
# Count and accumulate the results.
my $count = 0;
while (my $r = $s->next_result) {
    push @papers, { id => $r->pmid, title => $r->title, abstract => $r->abstract };
    $count++;
    print "$count papers stored.\n" if ($count % 100) == 0;
}
# Summarize the number of results.
print "$count total results.\n";
$vars{result_count} = $count;
# We have now compiled the information we need for the output. Create the template engine.
print "Constructing web page in $outFile.\n";
my $templateEngine = Template->new(ABSOLUTE => 1);
# Allocate the result variable.
my $html;
# Create the web page.
my $template = $opt->template;
print "Template file is $template.\n";
$templateEngine->process($opt->template, \%vars, \$html) || die $templateEngine->error();
# Output the web page.  For now, we kill wide characters without worrying about it.
$html =~ s/[^[:ascii:]]+//g;
print $oh $html;


## TODO process the input to produce the output.