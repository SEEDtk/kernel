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


=head1 Generate Synthetic Metagenome Data

    syn_mg_samples_pipeline [ options ] dataDir

This script generates synthetic reads for L<bins_create.pl> tests. The synthetic reads are generated according to instructions
present in an input file.

=head2 Parameters

There is one positional parameter-- the name of the data directory. This directory must contain an input file named C<input>.
The output is placed in subdirectories called C<Genomes> and C<Samples>.

The command-line options are those found in L<Shrub/script_options>.

=head2 Input File Format

The input file is a text file. Lines beginning with a pound sign (C<#>) are treated as comments.

The file contains three sections. Each section is preceded by a header that begins in the first column and consists of the
section name followed by a colon (C<:>).

=over 4

=item Genomes

This section should contain one line for each source genome. Each line contains a tab, the genome ID, one or more whitespace
characters, and the genome name.

=item Samples

This section should contain one line for each sample. Each line contains a tab, a sample ID, a hyphen (surrounded by
optional whitespace) and then one or more coverage specifications, space-delimited. Each coverage specification
consists of a genome ID, a comma, and a level of coverage.

=item Reads

This section should contain one line for each sample. Each line contains a tab, a sample ID, a hyphen (surrounded by
optional whitespace) and then one or more sampling specifications, space-delimited. Each sampling specification consists
of a genome ID, a comma, a read length, a second comma, and a percent mutation rate.

=back

So, for example, the following file generates data from E Coli and Strep pyogenes.

  Genomes:
        83333.1  Escherichia coli K12
        160490.1 Streptococcus pyogenes M1 GAS

  Samples:
        1 - 83333.1,20 160490.1,50
        2 - 83333.1,5  160490.1,100

  Reads:
        1 - 83333.1,500,5  160490.1,500,5
        2 - 83333.1,500,5  160490.1,500,5

The first sample has E coli at coverage 20 and Strep pyogenes at coverage 50. The second sample has E coli at coverage 5
and Strep pyogenes at coverage 100. When generating the reads, each read will have a length of 500 and a 5% mutation rate.

=head2 Output Files

There are several working directories and files produced. The primary output files are called C<all.reads>, underneath the
C<Samples> subdirectory of the data directory specified on the command line. Each sample's read file will be in a subdirectory
of C<Samples> with the same name as the sample ID number.

=cut

use strict;
use Data::Dumper;
use Carp;
use SeedUtils;
use gjoseqlib;
use Shrub;
use ScriptUtils;
use File::Copy::Recursive;


my $opt = ScriptUtils::Opts("dataDir", Shrub::script_options(),
            ['post', 'post-process the reads into assemblies']);
my $dataD = shift @ARGV;
my $shrub = Shrub->new_for_script($opt);
if (! -s "$dataD/input")  { die "usage: syn_mg_samples_pipeline DataD" }

my $input   = &read_input("$dataD/input");
print "Copying genomes.\n";
&generate_genomes($dataD,$input);

my $samples = $input->{Samples};
foreach my $sampleS (sort { $a->[0] <=> $b->[0] } @$samples)
{
    my $sample = $sampleS->[0];
    print "Generating sample $sample.\n";
    &generate_data_for_sample($sampleS,$input,$dataD);
    # Combine all reads into all.reads.
    print "Combining reads for $sample.\n";
    open(my $oh, ">$dataD/Samples/$sample/all.reads") || die "Could not open all.reads for $sample: $!";
    my $sampleDir = "$dataD/Samples/$sample/Reads";
    openDir(my $dh, $sampleDir) || die "Could not open reads directory for $sample: $!";
    my @readFiles = grep { $_ =~ /^\d+\.\d+$/ } readdir $dh;
    for my $readFile (@readFiles) {
        open(my $ih, "<$sampleDir/$readFile") || die "Could not open read file $readFile in $sampleDir: $!";
        while (! eof $ih) {
            my $line = <$ih>;
            print $oh, $line;
        }
    }
}

sub read_input {
    my($file) = @_;

    my $input = {};
    open(INPUT,"<") || die "could not open $file";
    $/ = "\n\n";
    while ($_ = <INPUT>)
    {
        $_ =~ s/^#.+?\n//gm;
        chomp;
        if ($_ =~ /^(\S+):\n\s*(\S.*\S)/s)
        {
            my $k = $1;
            my $x = $2;

            if ($k eq 'Genomes')
            {
                my @genomes = map { ($_ =~ /^\s*(\d+\.\d+)/) ? $1 : () } split(/\n/,$x);
                $input->{'Genomes'} = \@genomes;
            }
            elsif ($k eq "Samples")
            {
                my @samples = map { ($_ =~ /^\s*(\d+)\s*-\s*(\S.*\S)/) ? [$1,&split_parms($2)] : () }
                split(/\n/,$x);
                $input->{'Samples'} = \@samples;
            }
            elsif ($k eq "Mutations")   ### NOT IMPLEMENTED (AND PROBABLY SHOULD NOT BE)
            {
                my @mutations = map { ($_ =~ /^\s*(\d+)\s*-\s*(\S.*\S)/) ? [$1,&split_parms($2)] : () }
                split(/\n/,$x);
                $input->{'Mutations'} = \@mutations;
            }
            elsif ($k eq "Reads")
            {
                my @reads = map { ($_ =~ /^\s*(\d+)\s*-\s*(\S.*\S)/) ? [$1,&split_parms($2)] : () }
                split(/\n/,$x);
                $input->{'Reads'} = \@reads;
            }
            elsif ($k eq "Breaks") ### NOT IMPLEMENTED (AND PROBABLY SHOULD NOT BE)
            {
                my @breaks = map { ($_ =~ /^\s*(\d+)\s*-\s*(\S.*\S)/) ? [$1,&split_parms($2)] : () }
                split(/\n/,$x);
                $input->{'Breaks'} = \@breaks;
            }
        }
    }
    return $input;
}

sub split_parms {
    my($p) = @_;
    my @parsed = map { [split(/\,/,$_)] }
    split(/\s+/,$p);
    return \@parsed;
}

sub generate_data_for_sample {
    my($sampleS,$input,$dataD) = @_;

    my($sample,$encoded) = @$sampleS;
    if (! -d "$dataD/Samples/$sample")
    {
        if (! -d "$dataD/Samples")
        {
            mkdir("$dataD/Samples",0777) || die "could not make $dataD/Samples";
        }
        mkdir("$dataD/Samples/$sample",0777) || die "could not make $dataD/Samples/$sample";
        print "Generating close genomes for $sample.\n";
        &generate_close_genomes($dataD,$sample,$input);
        print "Generating reads genomes for $sample.\n";
        &generate_reads($dataD,$sample,$input);
    }
}

sub generate_genomes {
    my($dataD,$input) = @_;

    mkdir("$dataD/Genomes",0777);
    my @genomes = sort { $a <=> $b } @{$input->{Genomes}};
    foreach my $g (@genomes)
    {
        if (mkdir("$dataD/Genomes/$g",0777))
        {
            &get_contigs($g,"$dataD/Genomes/$g/contigs");
            &get_peg_locs($g,"$dataD/Genomes/$g/peg.locs");
        }
    }
}

sub generate_close_genomes {
    my($dataD,$sample,$input) = @_;
    my $breaks    = &breaks($sample,$input);
    my $mutations = &mutations($sample,$input);
    mkdir("$dataD/Samples/$sample/CloseGenomes",0777);
    my @genomes = sort { $a <=> $b } @{$input->{Genomes}};
    foreach my $g (@genomes)
    {
        &get_contigs($g,"$dataD/Samples/$sample/CloseGenomes/$g");
    }
}

sub generate_reads {
    my($dataD,$sample,$input) = @_;
    my $reads      = &reads($sample,$input);
    my $readsD     = "$dataD/Samples/$sample/Reads";
    mkdir($readsD,0777);
    my @genomes = sort { $a <=> $b } @{$input->{Genomes}};
    foreach my $g (@genomes)
    {
        my $coverage = &coverage($sample,$g,$input);
        my @close_genome_contigs = &gjoseqlib::read_fasta("$dataD/Samples/$sample/CloseGenomes/$g");
        open(READS,">$readsD/$g") || die "could not open $readsD/$g";
        my @tmp           = grep { $_->[0] eq $g } @$reads;
        my $len           = $tmp[0]->[1];
        my $mutation_rate = $tmp[0]->[2];
        &make_mutated_reads($sample,\@close_genome_contigs,$len,$mutation_rate,$coverage,\*READS);
        close(READS);
    }
}

sub make_mutated_reads {
    my($sample,$contigs,$len,$mutation_rate,$coverage,$fh) = @_;
    foreach my $tuple (@$contigs)
    {
        my($contig_id,undef,$seq) = @$tuple;
        my $contig_len = length($seq);
        my $last_pos = $contig_len - $len;
        my $num_reads = int(($contig_len * $coverage) / $len);
        my $nxt = 1;
        &write_mutated_read($sample,$contig_id,\$seq,0,$len,$nxt++,$fh,$mutation_rate);
        &write_mutated_read($sample,$contig_id,\$seq,$last_pos,$len,$nxt++,$fh,$mutation_rate);
        for (my $i=0; ($i < ($num_reads-2)); $i++)
        {
            my $pos = int(rand() * ($last_pos-1));
            &write_mutated_read($sample,$contig_id,\$seq,$pos,$len,$nxt++,$fh,$mutation_rate);
        }
    }
}

sub write_mutated_read {
    my($sample,$contig_id,$seqP,$pos,$len,$idN,$fh,$mutation_rate) = @_;

    my $seq = substr($$seqP,$pos,$len);
    my $mutations = int($mutation_rate * $len / 100);
    for (my $i=0; ($i < $mutations); $i++)
    {
        my $p = int(rand() * $len);
        my $c = &mutated(substr($seq,$p,1));
        substr($seq,$p,1) = $c;
    }
    print $fh ">$sample:$contig_id:$idN\n",$seq,"\n";
}

# maintains GC ratio
sub mutated {
    my($c) = @_;

    if    (($c eq 'a') || ($c eq 'A'))
    {
        return 'T';
    }
    elsif (($c eq 'c') || ($c eq 'C'))
    {
        return 'G';
    }
    elsif (($c eq 'g') || ($c eq 'G'))
    {
        return 'C';
    }
    else
    {
        return 'A';
    }
}

sub length_of_contigs {
    my($contigs) = @_;

    my $tot = 0;
    foreach my $tuple (@$contigs)
    {
        $tot += length($tuple->[2]);
    }
    return $tot;
}

sub coverage {
    my($sample,$genome,$input) = @_;

    my @tmp = grep { $_->[0] == $sample } @{$input->{Samples}};
    @tmp    = grep { $_->[0] eq $genome } @{$tmp[0]->[1]};
    return $tmp[0]->[1];
}

sub reads {
    my($sample,$input) = @_;

    my @tmp = grep { $_->[0] == $sample } @{$input->{Reads}};
    return $tmp[0]->[1];
}

sub breaks {
    my($sample,$input) = @_;

    my @tmp = grep { $_->[0] == $sample } @{$input->{Breaks}};
    return (@tmp > 0) ? $tmp[0]->[1] : [];
}

sub mutations {
    my($sample,$input) = @_;

    my @tmp = grep { $_->[0] == $sample } @{$input->{Mutations}};
    return (@tmp > 0) ? $tmp[0]->[1] : [];
}

##################################  change for SEEDtk

sub get_contigs {
    my($genome,$file) = @_;
    my $fname = $shrub->genome_fasta($genome);
    File::Copy::Recursive($fname, $file);
}

sub get_peg_locs {
#    my($genome,$file) = @_;
#    open(my $oh, ">$file") || die "Could not open peg output file $file: $!";
#    my @pegs = $shrub->GetFlat('Feature', 'Feature(id) LIKE ?', ["fig|$genome.peg.%"], 'id');
#    for my $peg (@pegs) {
#        my $loc = $shrub->loc_of($peg);
#        print $oh join("\t", $peg, $loc->String) . "\n";
#    }
}



