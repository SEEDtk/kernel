#!/usr/bin/env /vol/ross/FIGdisk/bin/run_perl

use strict;
use Data::Dumper;
use Carp;
use SeedUtils;
use gjoseqlib;
use Shrub;
use ScriptUtils;

my $opt = ScriptUtils::Opts("dataDir", Shrub::script_options());
my $dataD = shift @ARGV;
my $shrub = Shrub->new_for_script($opt);
if (! -s "$dataD/input")  { die "usage: syn_mg_samples_pipeline DataD" }

my $input   = &read_input("$dataD/input");
&generate_genomes($dataD,$input);

my $samples = $input->{Samples};
foreach my $sampleS (sort { $a->[0] <=> $b->[0] } @$samples)
{
    my $sample = $sampleS->[0];
    &generate_data_for_sample($sampleS,$input,$dataD);
    &SeedUtils::run("cat $dataD/Samples/$sample/Reads/* > $dataD/Samples/$sample/all.reads");
}

sub read_input {
    my($file) = @_;

    my $input = {};
    open(INPUT,"grep -v '^#' $file |") || die "could not open $file";
    $/ = "\n\n";
    while ($_ = <INPUT>)
    {
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
            elsif ($k eq "Mutations")
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
            elsif ($k eq "Breaks")
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
        &generate_close_genomes($dataD,$sample,$input);
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
        &write_mutated_read($sample,$contig_id,\$seq,0,$len,$nxt++,$fh);
        &write_mutated_read($sample,$contig_id,\$seq,$last_pos,$len,$nxt++,$fh);
        for (my $i=0; ($i < ($num_reads-2)); $i++)
        {
            my $pos = int(rand() * ($last_pos-1));
            &write_mutated_read($sample,$contig_id,\$seq,$pos,$len,$nxt++,$fh);
        }
    }
}

sub write_mutated_read {
    my($sample,$contig_id,$seqP,$pos,$len,$idN,$fh) = @_;

    print $fh ">$sample:$contig_id:$idN\n",substr($$seqP,$pos,$len),"\n";
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
    &SeedUtils::run("cp $fname $file");
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



