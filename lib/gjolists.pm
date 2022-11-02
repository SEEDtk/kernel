package gjolists;

#  Invoke with:
#
#     use gjolists;
#
#  List comparisons:
#
#  @common = common_prefix( \@list1, \@list2 )
#  @common = common_prefix_n( \@list1, \@list2, ... )
#  ( \@pref, \@suf1, \@suf2 ) = common_and_unique( \@list1, \@list2 )
#  ( \@suf1, \@suf2 )         = unique_suffixes( \@list1, \@list2 )
#
#  List properties:
#
#  @unique = unique_set( @list )     #  Reduce a list to a set
#  @dups   = duplicates( @list )
#
#  @random = random_order( @list )
#
#  Set algebra:
#
#  @A_or_B  = union( \@list1, \@list2, ... )
#  @A_and_B = intersection( \@list1, \@list2, ... )
#  @A_not_B = set_difference( \@list1, \@list2 )

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
        common_prefix
        common_prefix_n
        common_and_unique
        unique_suffixes
        
        unique_set
        duplicates
        random_order

        union
        intersection
        set_difference
        );

use strict;


#-----------------------------------------------------------------------------
#  Return the common prefix of two lists:
#
#  @common = common_prefix( \@list1, \@list2 )
#-----------------------------------------------------------------------------
sub common_prefix {
    my ($l1, $l2) = @_;
    ref($l1) eq "ARRAY" || die "common_prefix: arg 1 is not an array ref\n";
    ref($l2) eq "ARRAY" || die "common_prefix: arg 2 is not an array ref\n";
    my $i = 0;
    my $l1_i;
    while ( defined( $l1_i = $l1->[$i] ) && $l1_i eq $l2->[$i] ) { $i++ }

    return @$l1[ 0 .. ($i-1) ];  # perl handles negative range
}


#-----------------------------------------------------------------------------
#  Return the common prefix of two or more lists:
#
#  @common = common_prefix_n( \@list1, \@list2, ... )
#-----------------------------------------------------------------------------
sub common_prefix_n {
    my $n = @_;
    $n > 1 || die "common_prefix: requires 2 or more arguments\n";
    for (my $j = 1; $j <= $n; $j++) {
        ref($_[$j-1]) eq "ARRAY" || die "common_prefix_n: arg $j is not an array ref\n";
    }

    my $l0 = $_[0];
    my $l0_i;
    my $i;
    for ( $i = 0; defined( $l0_i = $l0->[$i] ); $i++ ) {
        for ( my $j = 1; $j < $n; $j++ ) {
            $l0_i eq $_[$j]->[$i] || ( return @$l0[0 .. ($i-1)] )
        }
    }

    return @$l0[ 0 .. ($i-1) ];  # perl handles negative range
}


#-----------------------------------------------------------------------------
#  Return the common prefix and unique suffixes of each of two lists:
#
#  ( \@prefix, \@suffix1, \@suffix2 ) = common_and_unique( \@list1, \@list2 )
#-----------------------------------------------------------------------------
sub common_and_unique {
    my ($l1, $l2) = @_;
    ref($l1) eq "ARRAY" || die "common_prefix: arg 1 is not an array ref\n";
    ref($l2) eq "ARRAY" || die "common_prefix: arg 2 is not an array ref\n";
    my $i = 0;
    my $l1_i;
    while ( defined( $l1_i = $l1->[$i] ) && $l1_i eq $l2->[$i] ) { $i++ }

    my $len1 = @$l1;
    my $len2 = @$l2;
    return ( [ @$l1[ 0  .. $i-1    ] ]  # perl handles negative range
           , [ @$l1[ $i .. $len1-1 ] ]
           , [ @$l2[ $i .. $len2-1 ] ]
           );
}


#-----------------------------------------------------------------------------
#  Return the unique suffixes of each of two lists:
#
#  ( \@suffix1, \@suffix2 ) = unique_suffixes( \@list1, \@list2 )
#-----------------------------------------------------------------------------
sub unique_suffixes {
    my ($l1, $l2) = @_;
    ref($l1) eq "ARRAY" || die "common_prefix: arg 1 is not an array ref\n";
    ref($l2) eq "ARRAY" || die "common_prefix: arg 2 is not an array ref\n";
    my $i = 0;
    my $l1_i;
    while ( defined( $l1_i = $l1->[$i] ) && $l1_i eq $l2->[$i] ) { $i++ }

    my $len1 = @$l1;
    my $len2 = @$l2;
    return ( [ @$l1[ $i .. $len1-1 ] ]  # perl handles negative range
           , [ @$l2[ $i .. $len2-1 ] ]
           );
}


#-----------------------------------------------------------------------------
#  Reduce a list to its unique elements (stable in order):
#
#  @unique = unique_set( @list )
#-----------------------------------------------------------------------------
sub unique_set {
    my %cnt = ();
    map { ( $cnt{$_} = $cnt{$_} ? $cnt{$_}+1 : 1 ) == 1 ? $_ : () } @_;
}


#-------------------------------------------------------------------------------
#  List of values duplicated in a list (stable in order by second occurance):
#
#  @dups = duplicates( @list )
#-------------------------------------------------------------------------------
sub duplicates
{
    my %cnt = ();
    grep { ++$cnt{$_} == 2 } @_;
}


#-------------------------------------------------------------------------------
#  Randomize the order of a list:
#
#  @random = random_order( @list )
#-------------------------------------------------------------------------------
sub random_order {
    my ( $i, $j );
    for ( $i = @_ - 1; $i > 0; $i-- ) {
        $j = int( ($i+1) * rand() );
        ( $_[$i], $_[$j] ) = ( $_[$j], $_[$i] );
    }

   @_
}


#-----------------------------------------------------------------------------
#  Union of two or more sets (by reference):
#
#  @union = union( \@set1, \@set2, ... )
#-----------------------------------------------------------------------------
sub union
{
    my %cnt = ();
    grep { ++$cnt{$_} == 1 } map { @$_ } @_;
}


#-----------------------------------------------------------------------------
#  Intersection of two or more sets:
#
#  @intersection = intersection( \@set1, \@set2, ... )
#-----------------------------------------------------------------------------
sub intersection
{
    my $set = shift;
    my @intersection = @$set;

    foreach $set ( @_ )
    {
        my %set = map { ( $_ => 1 ) } @$set;
        @intersection = grep { exists $set{ $_ } } @intersection;
    }

    @intersection;
}


#-----------------------------------------------------------------------------
#  Elements in set 1, but not set 2:
#
#  @difference = set_difference( \@set1, \@set2 )
#-----------------------------------------------------------------------------
sub set_difference
{
    my ($set1, $set2) = @_;
    my %set2 = map { ( $_ => 1 ) } @$set2;
    grep { ! ( exists $set2{$_} ) } @$set1;
}


1;
