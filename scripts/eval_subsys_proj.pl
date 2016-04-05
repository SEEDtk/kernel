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
use ScriptUtils;
use Shrub;
use SeedUtils;
use Data::Dumper;
use gjoseqlib;

=head1 Test Projection of Subsystems

    eval_subsys_proj -s subsystem -e errors -cmd CmpToAnnotate > report

This script evaluates the quality of projection for a subsystem on the coreSEED.

=head2 Output

The script produces two outputs:

=over 4

=item error-file

a file of 4-tuples of the form [peg,func-in-coreSEED,protein-sequence,projected-function]

=item report

three lines written to STDOUT containing

    matched     Matched
    mismatched  MisMatched
    success     Fraction-mismatched

=item Command_to_Annotate

    The name of an annotation script that is run using

        cmd training-file test-file > projections [feature\tcoreSEEDFunction\tsequence\tProjectedFunction]
        

=back

=cut

my $opt = ScriptUtils::Opts('', Shrub::script_options(),
          Shrub::script_options(), ScriptUtils::ih_options(),
          [ 'subsystem|s=s', 'Name or ID of a Subsystem on coreSEED', { required => 1 } ],
          [ 'errors|e=s','Name of created file of errors',{ required => 1 }],
          [ 'kmer_size|k=s','kmer size',{ default => 8 }],
          [ 'all_seqs|A=s','Created file of all sequences',{ }],
          [ 'cv|c=f','training fraction in cross-validation',{ default => 0.8 }],
          [ 'annotation_cmd|a=s','AnnotationCmd',{ required => 1 }]
    );
my $ih         = ScriptUtils::IH($opt->input);
my $shrub = Shrub->new_for_script($opt);

my $cmd        = $opt->annotation_cmd;
my $allF       = $opt->all_seqs;
my $cv         = $opt->cv;

my $k          = $opt->kmer_size;
my $error_file = $opt->errors;
my $subsys     = $opt->subsystem;
if (! $subsys) { die "you need to specify a subsystem"}
my $name;

my @tmp        = `svc_echo $subsys | get_entities Subsystem --fields name`;
if (@tmp == 0)
{
    @tmp = `all_entities Subsystem | grep "$subsys"`;
    if (@tmp == 0)
    {
	die "Invalid subsystem: $subsys";
    }
}
chomp;
($subsys,$name) = split(/\t/,$tmp[0]);

# my @tuples = $shrub->GetAll("Subsystem2Role 
#                              Role 
#                              Role2Function 
#                              Function 
#                              Function2Feature 
#                              Feature 
#                              Feature2Genome 
#                              Genome AND
#                              Feature2Protein Protein",
# 			    "(Subsystem2Role(from-link) = ?) AND (Genome(core) = ?) AND (Function2Feature(security) = ?)",
# 			    [$subsys,1,2],
# 			    "Feature(id) Function(description) Protein(sequence)");

my @tuples = $shrub->GetAll("Subsystem2Role 
                              Role 
                              Role2Function 
                              Function 
                              Function2Feature Feature Feature2Genome Genome",
 			    "(Subsystem2Role(from-link) = ?) AND (Genome(core) = ?) AND (Function2Feature(security) = ?)",
			    [$subsys,1,2],
			    "Feature(id) Function(description)");

my %tuples = map { ($_->[0] => $_) } @tuples;  # now remove duplicates
@tuples    = map { $tuples{$_} } keys(%tuples);

my %feat_to_func = map { ($_->[0] => $_->[1]) } @tuples;

open(PEGS,">tmp.$$");
foreach my $tuple (@tuples)
{
    print PEGS join("\t",@$tuple),"\n";
}
close(PEGS);
&run("svc_fasta -c 1 --prot < tmp.$$ > tmp.$$.fasta");
my @seq = &gjoseqlib::read_fasta("tmp.$$.fasta");
unlink("tmp.$$","tmp.$$.fasta");

foreach $_ (@seq)
{
    $_->[1] = $feat_to_func{$_->[0]};
}

if ($allF) {
    &gjoseqlib::write_fasta($allF,\@seq);
}

my @shuffled = &randomize(\@seq);
my $sz_training = int($cv * @shuffled);
if ($sz_training == @shuffled) { $sz_training-- }
my $sz_test     = @shuffled - $sz_training;
print STDERR "size training = $sz_training\n";
print STDERR "size test      = $sz_test\n";

my $good = 0;
my $bad  = 0;

open(ERRS,">$error_file") || die "could not open $error_file";
my $i = 0;
while ($i < @shuffled)
{
    open(TRAINING,">tmp.training.$$") || die "could not open tmp.training.$$";
    open(TEST,">tmp.test.$$") || die "could not open tmp.test.$$";
    my $j = 0;
    while ($j < @shuffled)
    {
	if (($j < $i) || ($j >= ($i + $sz_test)))
	{
	    print TRAINING join("\t",@{$shuffled[$j]}),"\n";
	}
	else
	{
	    print TEST join("\t",@{$shuffled[$j]}),"\n";
	}
	$j++;
    }
    close(TRAINING);
    close(TEST);
    
    &run("$cmd tmp.training.$$ tmp.test.$$ > tmp.projections.$$ $k");
#   print &Dumper("tmp.test.$$","tmp.training.$$","tmp.projections.$$"); die "HERE";
    open(PROJECTIONS,"<tmp.projections.$$") || die "could not open tmp.projections.$$";
    while (defined($_ = <PROJECTIONS>))
    {
	chomp;
	my($peg,$func,$seq,$projected) = split(/\t/,$_);
	if ($func eq $projected)
	{
	    $good++;
	}
	else
	{
	    print ERRS $_,"\n";
	    $bad++;
	}
    }	
    close(PROJECTIONS);
    $i += $sz_test;
}
unlink("tmp.test.$$","tmp.training.$$"."tmp.prijections.$$");

close(ERRS);
print "matched\t$good\n";
print "mismatched\t$bad\n";
print "success\t",sprintf("%0.3f",$good/($good+$bad)),"\n";

use List::Util qw(shuffle);
sub randomize {
    my($list) = @_;
    my @random = &shuffle(@$list);
    return @random;
}

sub run {
    my($cmd) = @_;

#    print STDERR "running: $cmd\n";                                                                                   
    my $rc = system($cmd);
    if ($rc)
    {
        die "$rc: $cmd failed";
    }
}
