package AVL;

#=============================================================================
# gjoavllib.pm
#-----------------------------------------------------------------------------
#  Beware that the package name does not match the file name.
#
#  Methods for AVL trees:
#
#     $avl  = AVL->new()
#     $avl  = AVL->new( \&cmp_func )
#     $avl  = AVL->new_numeric()
#     $avl  = AVL->new_text()
#
#     $bool = $avl->add( $key )         # insert new key into tree
#     $tree = $avl->del( $key )         # delete a key from tree
#     $bool = $avl->exists( $key )      # true if the key is in the tree
#
#     $key  = $avl->first()             # first key in tree
#     $key  = $avl->last()              # last key in tree
#     $key  = $avl->pop_first()         # first key in tree, deleting it
#     $key  = $avl->pop_last()          # last key in tree, deleting
#     $key  = $avl->next_key( $value )  # first key > value
#     $key  = $avl->prev_key( $value )  # first key < value
#
#     $n    = $avl->size()              # number of keys in the tree
#     @keys = $avl->flatten()           # return ordered list of keys
#
#  Two internal functions that really should not be used:
#
#     $node = $avl->root()
#    \&func = $avl->cmp_func()
#
#=============================================================================

use strict;
use Data::Dumper;

#-----------------------------------------------------------------------------
#  Get a new AVL tree
#
#    $avl = AVL::new();                # defaults to text
#    $avl = AVL::new( \&cmp_func );    # user-supplied sort function
#
#    $avl = AVL::new_caseless();       # case-insensitive text
#    $avl = AVL::new_numeric();        # numeric
#    $avl = AVL::new_text();           # text (same as new)
#
#  For example, you could get a case insensitive text sort with:
#
#   $avl = AVL->new( sub { lc $_[0] cmp lc $_[1] } )
#
#-----------------------------------------------------------------------------
sub new
{
    shift if UNIVERSAL::isa( $_[0], __PACKAGE__ );

    my $cmpfunc = shift || sub { $_[0] cmp $_[1] };

    ref($cmpfunc) eq 'CODE' ? bless [ undef, $cmpfunc ], 'AVL' : undef;
}

sub new_caseless { AVL::new( sub { lc $_[0] cmp lc $_[1] } ) }
sub new_numeric  { AVL::new( sub {    $_[0] <=>    $_[1] } ) }
sub new_text     { AVL::new( sub {    $_[0] cmp    $_[1] } ) }


#-----------------------------------------------------------------------------
#  Insert new key into tree, returning success indication
#
#     $bool = $avl->add( $key )
#
#-----------------------------------------------------------------------------
sub add
{
    my ( $tree, $key ) = @_;
    return undef unless $tree and defined( $key );
    my ( $root, $done ) = add_1( $tree->root(), $key, $tree->cmp_func() );

    $done && $tree->set_root( $root ) ? 1 : 0;
}

#
#  ( $root_node, $added ) = add_1( $node, $key, $cmp_func )
#
sub add_1
{
    my ( $node, $key, $cmp_func ) = @_;
    return ( new_tip( $key ), 1 ) if ! $node;

    my $added = 0;
    my $dir = &$cmp_func( $key, key( $node ) );

    if ( $dir < 0 )
    {
        my $nl = l( $node );
        if ( ! $nl )
        {
            set_l( $node, new_tip( $key ) );
            update_h( $node );
            $added = 1;
        }
        else
        {
            my $n2;
            ( $n2, $added ) = add_1( $nl, $key, $cmp_func );
            $node = balance( $n2, r( $n2 ), $node );
        }
    }

    elsif ( $dir > 0 )
    {
        my $nr = r( $node );
        if ( ! $nr )
        {
            set_r( $node, new_tip( $key ) );
            update_h( $node );
            $added = 1;
        }
        else
        {
            my $n4;
            ( $n4, $added ) = add_1( $nr, $key, $cmp_func );
            $node = balance( $node, l( $n4 ), $n4 );
        }
    }
    #  Otherwise it is not new, so we do not need to do anything

    ( $node, $added );
}


#-----------------------------------------------------------------------------
#  Delete a key from tree, returning false if it is not found.
#
#     $bool = $avl->del( $key )
#
#-----------------------------------------------------------------------------

sub del
{
    my ( $tree, $key ) = @_;

    my ( $node, $found ) = del_1( $tree->root(), $key, $tree->cmp_func() );
    $tree->set_root( $node ) if $found;
    $found;
}

sub del_1
{
    my ( $node, $key, $cmp_func ) = @_;
    return ( $node, 0 ) if ! $node;

    my $found = 0;
    my $dir = &$cmp_func( $key, key( $node ) );

    if ( $dir < 0 )
    {
        my $nl;
        ( $nl, $found ) = del_1( l( $node ), $key, $cmp_func );
        if ( $found )
        {
            set_l( $node, $nl );
            my $n4 = r( $node );
            $node = balance( $node, l( $n4 ), $n4 );
        }
    }

    elsif ( $dir > 0 )
    {
        my $nr;
        ( $nr, $found ) = del_1( r( $node ), $key, $cmp_func );
        if ( $found )
        {
            set_r( $node, $nr );
            my $n2 = l($node);
            $node  = balance( $n2, r( $n2 ), $node );
        }
    }

    else                                              #  Found it
    {
        $node = join_nodes( l( $node ), r( $node ) );
        $found = 1;
    }

    ( $node, $found );
}


#-----------------------------------------------------------------------------
#  Return 1 if key is in tree, otherwise 0.
#-----------------------------------------------------------------------------
sub exists
{
    my ( $tree, $key ) = @_;
    my $node = key_node( $tree->root(), $key, $tree->cmp_func() );
    $node ? 1 : 0;
}


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  Return tree node with given key, or undef if not found.
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub key_node
{
    my ( $node, $key, $cmp_func ) = @_;
    return undef if ! $node;

    my $dir = &$cmp_func( $key, key( $node ) );

    $dir < 0 ? key_node( l( $node ), $key, $cmp_func ) :
    $dir > 0 ? key_node( r( $node ), $key, $cmp_func ) :
               $node;
}


#-----------------------------------------------------------------------------
#  Get the key (and value) of the leftmost tip
#
#      $key = $avl->first()
#      $key = $avl->last()
#
#  Remove first or last tip from tree and return its key and value.
#
#      $key = $avl->pop_first()
#      $key = $avl->pop_last()
#
#-----------------------------------------------------------------------------
sub first
{
    my $tree = shift;
    key( first_node( $tree->root() ) );
}


sub last
{
    my $tree = shift;
    key( last_node( $tree->root() ) );
}


sub pop_first
{
    my $tree = shift;
    my ( $node, $root ) = pop_first_node( $tree->root() );
    $tree->set_root( $root );
    key( $node );
}


sub pop_last
{
    my $tree = shift;
    key( pop_last_node( $tree->root() ) );
    my ( $node, $root ) = pop_last_node( $tree->root() );
    $tree->set_root( $root );
    key( $node );
}


#-----------------------------------------------------------------------------
#  Return the next key > $key in the tree, or undef if there is none.
#  Return the prev key < $key in the tree, or undef if there is none.
#
#     $nextkey  = next_key( $key, $tree )
#     $prev_key = prev_key( $key, $tree )
#
#-----------------------------------------------------------------------------
sub next_key
{
    my ( $tree, $key ) = @_;
    next_key_1( $tree->root(), $key, $tree->cmp_func(), undef )
}

sub next_key_1
{
    my ( $node, $key, $cmp_func, $nextkey ) = @_;
    return $nextkey if ! $node;

    my $nodekey = key( $node );
    my $dir = &$cmp_func( $key, $nodekey );
    $dir < 0 ? next_key_1( l( $node ), $key, $cmp_func, $nodekey )
             : next_key_1( r( $node ), $key, $cmp_func, $nextkey );
}


sub prev_key
{
    my ( $tree, $key ) = @_;
    prev_key_1( $tree->root(), $key, $tree->cmp_func(), undef )
}

sub prev_key_1
{
    my ( $node, $key, $cmp_func, $prevkey ) = @_;
    return $prevkey if ! $node;

    my $nodekey = key( $node );
    my $dir = &$cmp_func( $key, $nodekey );
    $dir > 0 ? prev_key_1( r( $node ), $key, $cmp_func, $nodekey )
             : prev_key_1( l( $node ), $key, $cmp_func, $prevkey );
}


#-----------------------------------------------------------------------------
#  Other basic properties:
#
#     $n    = $avl->size()
#     @keys = $avl->flatten()
#
#-----------------------------------------------------------------------------
sub size
{
    my $tree = shift;
    size_1( $tree->root() );
}


sub size_1 { $_[0] ? 1 + size_1( l($_[0]) ) + size_1( r($_[0]) ) : 0 }


sub flatten
{
    my $tree = shift;
    flatten_1( $tree->root() )
}


sub flatten_1
{
    ( local $_ = shift ) ? ( flatten_1( l($_) ), key($_), flatten_1( r($_) ) )
                         : ();
}


#=============================================================================
#  Internally used tree functions:
#=============================================================================
#
#   $tree = [ $rootnode, $compfunc ]
#
#   $rootnode = $avl->root();
#   $cmpfunc  = $avl->cmp_func();
#
#   $node     = $avl->set_root( $node );
#
#-----------------------------------------------------------------------------

sub root     { $_[0] && ref($_[0]) eq "AVL" ? $_[0]->[0] : undef }
sub cmp_func { $_[0] && ref($_[0]) eq "AVL" ? $_[0]->[1] : undef }

sub set_root { $_[0] && ref($_[0]) eq "AVL" ? ( $_[0]->[0] = $_[1] ) : undef }


#=============================================================================
#  Internally used node functions:
#=============================================================================
#
#     $node = [ $lref, $key, $rref, $h ]
#
#  Extract key, left, right, or height
#
#     $key    = key( $node )
#     $l_node = l( $node )
#     $r_node = r( $node )
#     $height = h( $node )
#
#-----------------------------------------------------------------------------

sub l       { $_[0] ? $_[0]->[0] : undef }
sub key     { $_[0] ? $_[0]->[1] : undef }
sub r       { $_[0] ? $_[0]->[2] : undef }
sub h       { $_[0] ? $_[0]->[3] : 0 }

sub set_l   { $_[0] ? ( $_[0]->[0] = $_[1] ) : undef }
sub set_key { $_[0] ? ( $_[0]->[1] = $_[1] ) : undef }
sub set_r   { $_[0] ? ( $_[0]->[2] = $_[1] ) : undef }
sub set_h   { $_[0] ? ( $_[0]->[3] = $_[1] ) : undef }


#-----------------------------------------------------------------------------
#  Is this a valid node?
#
#    $bool = node( $node )
#
#-----------------------------------------------------------------------------
sub node
{
    ( local $_ = shift )
        and ref($_) eq "ARRAY" && defined($_->[1]) && defined($_->[3])
}


#-----------------------------------------------------------------------------
#  Return first (leftmost) node, or last (rightmost) node.
#
#    $first = first_node( $node )
#    $last  = last_node(  $node )
#
#-----------------------------------------------------------------------------
sub first_node
{
    local $_ = shift or return undef;
    my $l = l( $_ );
    $l ? first_node( $l ) : $_;
}


sub last_node
{
    local $_ = shift or return undef;
    my $r = r( $_ );
    $r ? last_node( $r ) : $_;
}


#-----------------------------------------------------------------------------
#  Pop and return first (leftmost) node, or last (rightmost) node.
#
#      $first          = pop_first_node( $node, $parent )
#    ( $first, $root ) = pop_first_node( $node, $parent )
#
#      $last           = pop_last_node(  $node, $parent )
#    ( $last,  $root ) = pop_last_node(  $node, $parent )
#
#-----------------------------------------------------------------------------
sub pop_first_node
{
    my ( $node, $parent ) = @_;
    return undef if ! $node;

    my ( $rest, $l_node );
    my $nl = l( $node );
    if ( ! $nl )                 # This is the leftmost node
    {
        #  Move the right subtree off the node, replacing node on the parent.
        #  If there is no parent, the caller must check r($node).
        #
        #       p                    p
        #      /                    /
        #  l_node   ==>   l_node + nr
        #      \
        #      nr
        #

        $rest = r( $node );
        set_r( $node, undef );
        $l_node = $node;
    }

    else
    {
        #
        #  Continue descent into left subtree. Upon returning, we have:
        #
        #                 p
        #                /
        #   l_node +  node
        #             /   \
        #            /     \
        #          nl       nr
        #                  /
        #                 n3
        #
        #  Because l_node has been removed from nl, we might need to rebalance.
        #
        $l_node = pop_first_node( $nl, $node )
            or die "pop_first: bad tree\n";

        my $nr = r( $node );
        my $n3 = l( $nr );
        $rest  = balance( $node, $n3, $nr );
    }

    set_l( $parent, $rest ) if $parent;

    wantarray ? ( $l_node, $rest ) : $l_node;
}


sub pop_last_node
{
    my ( $node, $parent ) = @_;
    return undef if ! $node;

    my ( $rest, $r_node );
    my $nr = r( $node );
    if ( ! $nr )                   # This is rightmost node
    {
        #  Move the left subtree off the node, transplanting it to the parent.
        #  If there is no parent, the caller must check r($node) to get the
        #  new root of the tree.
        #
        #   p         p
        #    \         \
        #     n  ==>    nl + n=rtip
        #    /
        #  nl
        #

        $rest = l( $node );
        set_l( $node, undef );
        $r_node = $node;
    }
    else
    {
        #
        #  Continue descent into right subtree. Upon returning, we have:
        #
        #                 p
        #                /
        #              node  + r_node
        #             /    \
        #            /      \
        #          nl       nr
        #            \ 
        #            n3
        #
        #  Because r_node has been removed from nr, we might need to rebalance.
        #
        $r_node = pop_last_node( $nr, $node )
            or die "pop_last: bad tree\n";

        my $nl = l( $node );
        my $n3 = r( $nl );
        $rest  = balance( $nl, $n3, $node );
    }

    set_r( $parent, $rest ) if $parent;

    wantarray ? ( $r_node, $rest ) : $r_node;
}


#-----------------------------------------------------------------------------
#  Join two subtrees for which common parent has been deleted
#-----------------------------------------------------------------------------
#
#                     /  \
#                   nl    nr
#
#-----------------------------------------------------------------------------
sub join_nodes
{
    my ( $nl, $nr ) = @_;

    return $nr if ! $nl;      # Correctly handles n3 = undef
    return $nl if ! $nr;

    my $node;

    #  If the left subtree is higher, pull its rightmost node to serve as the
    #  new root.
    if ( h( $nl ) >= h( $nr ) )
    {
        my $rest;
        ( $node, $rest ) = pop_last_node( $nl, undef );
        set_l( $node, $rest );
        set_r( $node, $nr );
    }

    #  Otherwise, pull the leftmost node for the right subtree.
    else
    {
        my $rest;
        ( $node, $rest ) = pop_first_node( $nr, undef );
        set_l( $node, $nl );
        set_r( $node, $rest );
    }

    update_h( $node );
    $node;
}


#-----------------------------------------------------------------------------
# root subtrees to maintain balance
#-----------------------------------------------------------------------------
#
#                      n2    n4
#                     /  .  .  \
#                   n1    n3    n5
#                        /  \
#                      n3l  n3r
#
#                            $h1 >= $h3       $h5 >= $h3
#   ! $n2        ! $n4       $h1 >= $h5       $h5 >= $h1         otherwise
#  --------     --------     -----------      -----------      --------------
#     n4           n2           n2                  n4               n3
#    /  \         /  \         /  \                /  \             /  \
#  n3    n5     n1    n3     n1    n4            n2    n5         n2    n4
#                                 /  \          /  \             / \    / \
#                               n3    n5      n1    n3         n1 n3l  n3r n5
#
#-----------------------------------------------------------------------------
sub balance
{
    my ( $n2, $n3, $n4 ) = @_;

    if ( ! $n2 )
    {
        if ( ! $n4 ) { return $n3 }
        set_l( $n4, $n3 );
        update_h( $n4 );
        return $n4;
    }

    if ( ! $n4 )
    {
        set_r( $n2, $n3 );
        update_h( $n2 );
        return $n2;
    }

    my $h1 = h( l( $n2 ) );
    my $h3 = h(    $n3   );
    my $h5 = h( r( $n4 ) );

    if ( $h1 >= $h3 && $h1 >= $h5 )
    {
        set_r( $n2, $n4 );
        set_l( $n4, $n3 );
        update_h( $n4 );
        update_h( $n2 );
        return $n2;
    }

    if ( $h5 >= $h3 && $h5 >= $h1 )
    {
        set_r( $n2, $n3 );
        set_l( $n4, $n2 );
        update_h( $n2 );
        update_h( $n4 );
        return $n4;
    }

    else
    {
        my $n3l = l( $n3 );
        my $n3r = r( $n3 );
        set_r( $n2, $n3l );
        set_l( $n3, $n2 );
        set_r( $n3, $n4 );
        set_l( $n4, $n3r );
        update_h( $n2 );
        update_h( $n4 );
        update_h( $n3 );
        return $n3;
    }
}


#-----------------------------------------------------------------------------
#  Make a new tip node
#
#    $node = new_tip( $key )
#
#-----------------------------------------------------------------------------
sub new_tip { defined($_[0]) ? [ undef, $_[0], undef, 1 ] : undef }


#-----------------------------------------------------------------------------
#  Update the height of a node after an edit of the tree
#
#    update_h( $node )
#
#-----------------------------------------------------------------------------
sub update_h
{
    local $_ = shift or return undef;
    set_h( $_, max2( h( l($_) ), h( r($_) ) ) + 1 );
}

sub max2 { $_[0] >= $_[1] ? $_[0] : $_[1] }

1;
