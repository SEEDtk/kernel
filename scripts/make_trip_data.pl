use strict;
use Data::Dumper;
use SeedUtils;
use ScriptUtils;
use Shrub;
use Triples;

# usage: make_triples_data -d TripDir

=pod

=head1 Compute Data to Support Conserved Triple

      make_trip_data -d TriplesDirectory

Triples are sets of three functional roles from subsystems 
that tend to occur close to one another on the chromsomes (contigs).

First we identify the functional roless from subsystems.  These are the
only roles allowed in triples.

We go through all of the genomes making lists ofgenes implementing these roles, 
sorting them by location, and then computing the sets of triples.

NOTE: Two genes implemeting these roles are considered contiguous
no matter what distance separates them, and interveing genes that are not
in subsystems are just ignored.

This program builds a directory of related data files.

=head2 Output

The directory must be specified as an input argument (-d Dir), and it 
will be created, if it does not already exist.  The following files are
constructed in the directory

=over 4

=item sorted.all.triples

This is a file containig all triples from all genomes.

=item filtered triples

This is the condensed set of triples, which includes the number of occurrences
for each triple.  Each line contains 

    [number occurrences,role1, role2, role3, peg1, peg2, peg3].

Only trples that occur more than min_hits (default 20) times are kept.

=item clusters

This is a file of clustered roles (two roles are considered connected
if they occur in a triple that passed the filtering, and clustering is based
on transitive closure (two connected roles will occur in the saame cluster).
A cluster is a tab-separated list of roles on a single line.

=item clusters.rel

Here clusters are given in a 3-column format 

    [clusterID,roleID,role description]

=back

=head2 parameters

There are no positional parameters. The standard input is presumed to contain a list of
genomes that are considered solid examples of the subsystem. The command-line options are
those found in L<Shrub/script_options> and L<ScriptUtils/ih_options> plus the following.

=over 4

=item -d Directory

The directory that gets build containing the emerging data relating to triples.

=item -m MinOcc

The minimum number of times a triple must exist for it to make the
set of filtered.triples.

=back

=cut


my $opt = ScriptUtils::Opts(
    '',
    Shrub::script_options(),
    ['tripdata|d=s','Data Directory for Triples',{ required => 1 }],
    ['min|m=s','Minimum number of hits',{ default => 20 }]
);
my $dataD = $opt->tripdata;
my $min_hits = $opt->min;
mkdir($dataD,0777);

my $shrub = Shrub->new_for_script($opt);
my @tuples = $shrub->GetAll("Genome", "", [], "Genome(id)");
my @genomes = map { $_->[0] } @tuples;


my $roles = &Triples::subsystem_roles($shrub);
my $counts = {};
if (! -s "$dataD/sorted.all.triples")
{
    open(TRPS,"| sort -T . > $dataD/sorted.all.triples")
	|| die "could not create sorted triples";

    foreach my $genome (sort { $a <=> $b } @genomes)
    {
	print STDERR "$genome=$genome\n";
	&Triples::triples_for_genome($roles,$genome,$shrub,$counts,\*TRPS);
    }
    close(TRPS);
}

if (! -s "$dataD/filtered.triples")
{
    open(ALL,"<$dataD/sorted.all.triples")     || die "could not open $dataD/sorted.all.triples";
    open(FILTERED,"| sort -k 1 -n -r -T . > $dataD/filtered.triples")  || die "could not open $dataD/filtered.triples";

    while (defined($_ = <ALL>))
    {
	chomp;
	my($r1,$p1,$r2,$p2,$r3,$p3) = split(/\t/,$_);
	my $hits = $counts->{join("\t",($r1,$r2,$r3))};
	if ($hits && ($hits >= $min_hits))
	{
	    print FILTERED join("\t",($hits,$r1,$r2,$r3)),"\t",join("\t",($p1,$p2,$p3)),"\n";
	}
    }
    close(ALL);
    close(FILTERED);
}

if (! -s "$dataD/clusters")
{
    open(TRP,"<$dataD/filtered.triples")
	|| die "could not open $dataD/filtered.triples";
    open(CLUST,"| cluster_objects > $dataD/clusters")
	|| die "could not cluster objects";
    while (defined($_ = <TRP>))
    {
	if ($_ =~ /^\d+\t(\S+)\t(\S+)\t(\S+)/)
	{
	    if ($1 ne $2) { print CLUST "$1\t$2\n" }
	    if ($1 ne $3) { print CLUST "$1\t$3\n" }
	    if ($2 ne $3) { print CLUST "$2\t$3\n" }
	}
    }
    close(TRP);
    close(CLUST);

    my $role_desc = &Triples::subsystem_roles($shrub);
    open(TOREL,"tabs2rel < $dataD/clusters |")
	|| die "could not convert to rel format";
    open(REL,">$dataD/clusters.rel") || die "could not open $dataD/clusters.rel";
    while (defined($_ = <TOREL>))
    {
	if ($_ =~ /^(\S+)\t(\S.*\S)/)
	{
	    print REL join("\t",($1,$2,$role_desc->{$2})),"\n";
	}
    }
    close(TOREL);
    close(REL);
}
