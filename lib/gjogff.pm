package gjogff;

use strict;
use gjoseqlib;
use Data::Dumper;

#-------------------------------------------------------------------------------
#  Parse the descriptions of the features into key-value pairs
#
#     ( \%features_by_contig, \@dna ) = read_gff( \*FH, \%options );
#
#  Options:
#
#     contig_ids => \%ids
#     contig_ids => \@ids
#     dna        => \@dna
#
#  Features are:
#
#     [ $type, $loc, $id, $name, \%keys_values ]
#
#     Location is Sapling style contig_beg±length
#
#-------------------------------------------------------------------------------
sub read_gff
{
    my $opts = @_ && $_[ 0] && ref $_[ 0] eq 'HASH' ? shift :
               @_ && $_[-1] && ref $_[-1] eq 'HASH' ? pop   :
                                                      {};

    my $contig_ids = $opts->{ contig_ids } && ref $opts->{ contig_ids } eq 'HASH'  ? $opts->{ contig_ids }                                     :
                     $opts->{ contig_ids } && ref $opts->{ contig_ids } eq 'ARRAY' ? { map { $_ => 1 } @{ $opts->{ contig_ids } } }            :
                     $opts->{ contig_len } && ref $opts->{ contig_len } eq 'HASH'  ? $opts->{ contig_len }                                     :
                     $opts->{ contig_len } && ref $opts->{ contig_len } eq 'ARRAY' ? { map { $_ => 1 } @{ $opts->{ contig_len } } }            :
                     $opts->{ dna }        && ref $opts->{ dna }        eq 'ARRAY' ? { map { $_->[0] => length $_->[2] } @{ $opts->{ dna } } } :
                                                                                     undef;

    my ( $fh, $close ) = input_filehandle( $_[0] );

    my $dna;
    my %ftrs;
    local $_;

    my %no_contig;
    while ( defined( $_ = <$fh> ) )
    {
        if ( m/^##FASTA/i ) { $dna = 1; last }  #  GFF DNA introduction
        next if m/^track/i;                     #  JGI header
        next if m/^#/;                          #  GFF comment

        chomp;
        my ( $contig, $auth, $type, $l, $r, $scr, $dir, $fr, $rest ) = split /\t/;
        my ( $id, $name, $key_value ) = parse_description( $rest );
        push @{ $ftrs{ $id }->{ $type } },
             [ $id, $name, $contig, $l, $r, $dir, $fr, $key_value, $l+$r ];
    }

    my @dna;
    @dna = gjoseqlib::read_fasta( $fh ) if $dna;
    $contig_ids ||= { map { $_->[0] => 1 } @dna }  if @dna;

    close( $fh ) if $close;

    #  Merge segments of same feature.  Resulting descriptions are:
    #
    #     [ $type, $loc, $id, $name, $frame, \%keys_values ]
    #

    my $features_by_contig = process_gff_features( \%ftrs, $contig_ids );

    wantarray ? ( $features_by_contig, \@dna ) : $features_by_contig;
}


#-------------------------------------------------------------------------------
#  Parse the descriptions of the features into key-value pairs
#
#     ( $id, $name, \%key_values ) = parse_description( $gff_last_field );
#
#-------------------------------------------------------------------------------
sub parse_description
{
    #  Check for %xx escaped text:

    my @candidates = $_[0] =~ m/(\%..)/gi;
    my $escaped = @candidates && ( @candidates == grep { m/\%[0-9a-f][0-9a-f]/i } @candidates );

    #  Split into key=value pairs, remove surounding quotes, and unescape %..
    #  Some keys can have multiple values, so keep lists.

    my %key_values;
    foreach ( split / *; */, $_[0] )
    {
        my ( $key, $value ) = /^(\S+)[ =](.+)$/ ? ( $1 => $2 ) : ( 'name' => $_ );
        $value =~ s/""/"/g if $value =~ s/^"(.*)"$/$1/;
        $value =~ s/''/'/g if $value =~ s/^'(.*)'$/$1/;
        $value = unescape( $value ) if $escaped;
        push @{$key_values{$key}}, $value;
    }

    #  Look for an ID and name

    my $key;

    ( $key ) = grep { lc $_ eq 'name' } keys %key_values;
    ( $key ) = grep { m/name/i }        keys %key_values if ! $key;
    ( $key ) = grep { /I[dD]$/ }        keys %key_values if ! $key; 
    ( $key ) = grep { /id$/ }           keys %key_values if ! $key;
    my $name = $key ? $key_values{ $key }->[0] : $_[0];

    ( $key ) = grep { m/^ID$/i }       keys %key_values; 
    ( $key ) = grep { /I[dD]$/ }       keys %key_values if ! $key; 
    ( $key ) = grep { /id$/ }          keys %key_values if ! $key;
    ( $key ) = grep { lc eq 'parent' } keys %key_values if ! $key;
    my $id = $key ? $key_values{ $key }->[0] : $name;

    ( $id, $name, \%key_values );
}


#-------------------------------------------------------------------------------
#  Integrate and sort feature information:
#
#     \%features_by_contig = process_gff_features( \%features, \%contig_ids );
#
#  Features are:
#
#     [ $type, $loc, $id, $name, \%keys_values ]
#
#     Location is Sapling style contig_beg±length
#
#-------------------------------------------------------------------------------
sub process_gff_features
{
    my $ftrs       = shift || {};
    my $contig_ids = $_[0] && ref $_[0] eq 'HASH' ? shift : {};
    my $have_ids   = keys %$contig_ids;
    my $have_len   = ( values %$contig_ids )[0] > 1;

    my @ftrs;
    my %missing_contig;              #  Not currently reporting
    foreach my $id ( keys %$ftrs )
    {
        my $types = $ftrs->{ $id };
        foreach my $type ( keys %$types )
        {
            my $okay = 1;
            my $data = $types->{ $type };

            #              0     1       2       3       4     5       6         7          8 
            #  $datum = [ $id, $name, $contig, $left, $right, $dir, $frame, \%keys_values, $mid ]

            #
            #  Contig anlysis:
            #
            #  Number of contigs  contig_data   ids_matching      response
            #  -----------------------------------------------------------------
            #          1              yes            0         contig missing
            #          1              yes            1              okay
            #          1              no                            okay
            #         >1              yes            0         contig missing
            #         >1              yes            1           keep the 1
            #         >1              yes           >1         multiple contigs
            #         >1              no                       multiple contigs
            #  -----------------------------------------------------------------
            #

            my %contigs     = map { $_->[2] => 1 } @$data;
            my @contigs     = keys %contigs;
            my @have_contig = grep { $contig_ids->{ $_ } } @contigs;
            if ( @contigs == 1 )
            {
                $missing_contig{ $contigs[0] } = 1 if $have_ids && ! @have_contig;
            }
            elsif ( @have_contig == 1 )
            {
                #  Filter raw data to the contig we have:
                @$data = grep { $_->[2] eq $have_contig[0] } @$data;
            }
            elsif ( $have_ids && ! @have_contig )
            {
                $missing_contig{ $contigs[0] } = 1;
            }
            else
            {
                print STDERR "Feature '$id' of type '$type' has multiple contigs: ",
                              join( ', ', sort @contigs ), "\n";
                $okay = 0;
            }

            my %names = map { $_->[1] => 1 } @$data;
            if ( keys %names > 1 )
            {
                print STDERR "Feature '$id' of type '$type' has multiple names: ",
                              join( ', ', sort keys %names ), "\n";
#                $okay = 0;
            }

            my %dirs = map { $_->[5] => 1 } @$data;
            if ( keys %dirs > 1 )
            {
                print STDERR "Feature '$id' of type '$type' has multiple directions.\n";
                $okay = 0;
            }

            next if ! $okay;

            #  Grab shared information from the first exon:

            my ( $name, $cont, $dir, $frame ) = @{ $data->[0] }[ 1, 2, 5, 6 ];

            #  Sort my midpoint

            my $ndir = $dir =~ /^-/ ? -1 : +1;
            my @segs = sort { $ndir * ( $a->[8] <=> $b->[8] ) } @$data;

            #  Location in Sapling format
 
            my @locs = map { [ $cont,                              # contig
                               ( $ndir > 0 ? $_->[3] : $_->[4] ),  # start
                               $_->[5],                            # dir
                               $_->[4] - $_->[3] + 1               # length
                             ] }
                       @segs;
            my $loc = join( ',', map { "$_->[0]_$_->[1]$_->[2]$_->[3]" } @locs );

            my $mid = 0.5 * ( $segs[0]->[3] + $segs[-1]->[4] );

            #  Merge the keys and values for the merged gff lines:

            my %keys_values;
            my %n_seen;
            foreach my $keys_values0 ( map { $_->[7] } @$data )
            {
                foreach my $key ( keys %$keys_values0 )
                {
                    my $values = $keys_values0->{ $key };
                    foreach ( @$values )
                    {
                        push @{ $keys_values{ $key } }, $_  if ! $n_seen{ $key }->{ $_ }++;
                    }
                }
            }

            if ( $have_len && $contig_ids->{ $cont } )
            {
                if ( max( map { $_->[4] } @$data ) > $contig_ids->{ $cont } )
                {
                    print STDERR "Feature $type: $id $name extends beyond end of the contig DNA.\n";
                    next;
                }
            }
            
            #  The JGI files list frame 0 in all forward genes, so kill this.
            #  $keys_values{ codon_start } = [ $frame =~ /\d/ ? $frame+1 : 1 ];

            push @ftrs, [ [ $type, $loc, $id, $name, \%keys_values ], $cont, $mid ];
        }
    }

    my %features_by_contig;
    foreach ( sort { $a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] } @ftrs )
    {
        push @{ $features_by_contig{ $_->[1] } }, $_->[0];
    }

    if ( $have_ids && keys %missing_contig )
    {
        print STDERR "No data for contigs:\n";
        print STDERR map { "    $_\n" } sort keys %missing_contig;
    }

    \%features_by_contig;
}


sub max { my $max = shift; foreach ( @_ ) { $max = $_ if $_ > $max } $max }


#-------------------------------------------------------------------------------
#  Deal with %.. escaped text
#
#     $text = unescape( $text );
#
#-------------------------------------------------------------------------------
sub unescape
{
    join '', map { /\%(..)/ ? chr( hex( $1 ) ) : $_ }  # convert
             $_[0] =~ m/(\%[0-9a-f][0-9a-f]|.)/gi;     # split
}


#-------------------------------------------------------------------------------
#  Helper function for defining an input filehandle:
#
#     filehandle is passed through
#     string is taken as file name to be openend
#     undef or "" defaults to STDOUT
#
#      \*FH           = input_filehandle( $file );
#    ( \*FH, $close ) = input_filehandle( $file );
#
#-------------------------------------------------------------------------------
sub input_filehandle
{
    my $file = shift;

    #  Null string or undef

    if ( ! defined( $file ) || ( $file eq '' ) )
    {
        return wantarray ? ( \*STDIN, 0 ) : \*STDIN;
    }

    #  FILEHANDLE?

    if ( ref( $file ) eq "GLOB" )
    {
        return wantarray ? ( $file, 0 ) : $file;
    }

    #  File name

    if ( ! ref( $file ) )
    {
        -f $file or die "Could not find input file \"$file\"\n";
        my $fh;
        open( $fh, "<$file" ) || die "Could not open \"$file\" for input\n";
        return wantarray ? ( $fh, 1 ) : $fh;
    }

    return wantarray ? ( \*STDIN, undef ) : \*STDIN;
}


1;
