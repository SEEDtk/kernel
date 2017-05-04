package RepKmers;

use strict;
use Data::Dumper;

sub sim {
    my($s1,$s2,$K) = @_;

    my %in1;
    my $common = 0;
    my $i;
    for ($i=0; ($i < (length($s1)-$K)); $i++)
    {
        $in1{lc substr($s1,$i,$K)} = 1;
    }

    for ($i=0; ($i < (length($s2)-$K)); $i++)
    {
        if ($in1{lc substr($s2,$i,$K)})
        {
            $common++;
        }
    }
    return $common;
}

sub kmers_of_seq {
    my($K,$seq) = @_;

    my $kmers = {};
    my $i;
    for ($i=0; ($i < (length($seq)-$K)); $i++)
    {
        $kmers->{lc substr($seq,$i,$K)} = 1;
    }
    return $kmers;
}

sub in_common {
    my($K,$kmers,$seq) = @_;

    my $common = 0;
    my $i;
    for ($i=0; ($i < (length($seq)-$K)); $i++)
    {
        if ($kmers->{lc substr($seq,$i,$K)})
        {
            $common ++;
        }
    }
    return $common;
}

# returns a list of 2-tuples [CommonSc,Index]
sub closest_N_genomes {
    my($cached,$index,$n) = @_;

    my @closest;
    my $sims = $cached->{sims};
    my $closest = $sims->[$index];
    for (my $i=0; ($i < @$closest) && ($i < $n); $i++)
    {
        push(@closest,$closest->[$i]);
    }
    return @closest;
}

# returns a list of 2-tuples [CommonSc,Index]
sub closest_by_sc {
    my($cached,$index,$sc) = @_;

    my @closest;
    my $sims = $cached->{sims};
    my $closest = $sims->[$index];
    for (my $i=0; ($i < @$closest) && ($closest->[$i]->[0] >=  $sc); $i++)
    {
        push(@closest,$closest->[$i]);
    }
    return @closest;
}

sub rep_guts {
    my($cached,$max_sim,$todo) = @_;

    my $sims = $cached->{sims};
    my %seen;
    my $reps = [];
    while (defined(my $i = shift @$todo))
    {
        if (! $seen{$i})
        {
            push(@$reps,$i);
            $seen{$i} = 1;
#           print STDERR "added $i to reps\n";
            if (my $close = $sims->[$i])
            {
                if (defined($close))
                {
                    my $i2 = 0;
                    while (($i2 < @$close) && ($close->[$i2]->[0] >= $max_sim))
                    {
                        $seen{$close->[$i2]->[1]} = 1;
#			print STDERR "marked $close->[$i2]->[1]\n";
                        $i2++;
                    }
                }
            }
        }
    }
    return @$reps;
}

sub find_closest {
    my($cached,$best_id,$best_so_far,$max_sim,$kmersQ) = @_;

    my $index_to_g = $cached->{index_to_g};
    my $sims = $cached->{sims};
    my $seqs = $cached->{seqs};
    my %seen;
    while ($max_sim < 200)
    {
        my @close_tuples = &closest_by_sc($cached,$best_id,$max_sim);
        foreach my $tuple (@close_tuples)
        {
            my($sc,$id) = @$tuple;
            if (! $seen{$id})
            {
                $seen{$id} = 1;
            }
        }
        my $todo = [keys(%seen)];
        my @reps = &RepKmers::rep_guts($cached,$max_sim * 2,$todo);
        foreach my $id (@reps)
        {
            my $gid     = $index_to_g->{$id};
            my $seq = $seqs->{$gid};
            my $sc = &in_common(8,$kmersQ,$seq);
            if ($sc > $best_so_far)
            {
                $best_id = $id;
                $best_so_far = $sc;
            }
        }
        $max_sim = 2 * $max_sim;
    }
    return ($best_id,$best_so_far);
}

sub extract_kmers {
    my($seq,$K) = @_;

    my $triples = int($K/2);
    my $len = 3 * $triples;
    my @kmers;
    my $pos = 0;
    my $last = length($seq) - $len;

    while ($pos <= $last)
    {
        my $kmer = substr($seq, $pos, $len);
        $kmer =~ s/(..)./$1/g;
        if ($kmer =~ /^[acgt]+$/) {
            push @kmers, $kmer;
        }
        $pos++;
    }
    return \@kmers;
}

#
#  $seq is a sequence (a contig)
#  $pos is an index of where to start pulling a k-mer
#  $triples is the number of "2of3" characters we pull.
#
#  The returned string is composed of $trples 2of3 chunks
#  Thus, a $triples value of 8 would pull 16 characters.
#
sub extract_kmer {
    my($seqP,$pos,$triples) = @_;

    my(@chars);
    my $i;
    for ($i=0; ($i < $triples); $i++)
    {
        my $kmer = lc substr($$seqP,$pos+($i*3),2);
        if ($kmer =~ /^[acgt]*$/) {
            push(@chars, $kmer);
            #push(@chars,lc substr($$seqP,$pos+($i*3),2));
        } else {
            print STDERR "Funny at position $pos $i.\n";
        }
    }
    return join("",@chars);
}

1;
