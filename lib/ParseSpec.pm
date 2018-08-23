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


package ParseSpec;

    use strict;
    use warnings;
    use File::Slurp;
    use Data::Dumper;
    use CGI;
    use Regexp::Common;

=head1 Parse a kBase Spec File for Types

This object contains methods for parsing a kBase type specification and producing HTML output. It contains the following
members.

=over 4

=item typeHash

A reference to a hash mapping each type name to its a 2-tuple consistion of (0) its definition and (1) a comment. A type
definition is a list reference consisting of a category (C<mapping>, C<tuple>, C<list>, C<structure>, or another type), and
zero or more elements. Each element consists of a type, a comment, and an optional member name. The member name is only
present for structures.

=item comment

The text of the last comment encountered.

=back

=head2 Special Methods

=head3 new

    my $parser = ParseSpec->new();

Create a new, empty parser.

=cut

sub new {
    my ($class) = @_;
    my %typeHash = ('int' => ['basic integer type', 'type'], 'float' => ['floating-point numeric value', 'type'],
                    'string' => ['variable-length text string', 'type']);
    my $retVal = { typeHash => \%typeHash, comment => '' };
    bless $retVal, $class;
    return $retVal;
}


=head2 Public Methods

=head3 ToHtml

    my $html = $parser->ToHtml($specFile, $cgi);

Parse a kBase specification file and output the type definitions as Html. Each type definition will
be given a section in the table of contents, and references to other types in the file will be
linked. Comment sections will be associated with the subsequent type definition.

This method will make two passes over the data. The first pass will break the document into
sections based on the defined type. The second pass will turn each defined type into HTML.

=over 4

=item specFile

Name of the file containing the specification.

=item cgi

L<CGI> object for producing HTML output.

=item RETURN

Returns an HTML string.

=back

=cut

sub ToHtml {
    # Get the parameters.
    my ($self, $specFile) = @_;
    # Slurp in the spec file.
    my $spec = File::Slurp::read_file($specFile);
    # Extract the module name.
    my ($modComment, $modName, $modText) =
            $spec =~ /(?:\/\*\s*(.+?)\s*\*\/)?\s*module\s+(\S+)\s+{(.+)};/s;
    if (! $modName) {
        $modName = "Anonymous";
        $modText = $specFile;
    }
    # First split out the comments.
    my @csplits = split /(\/\*.+?\*\/)/s, $modText;
    # Now loop through the text, extracting tokens.
    my @tokens;
    for my $section (@csplits) {
        if ($section =~ /\/\*/) {
            push @tokens, $section;
        } else {
            push @tokens, grep { $_ } split /\s+|([<>{},;])/s, $section;
        }
    }
    # Now we need to loop through the tokens, creating type definitions.
    while (@tokens) {
        $self->ParseTokens(\@tokens)
    }
    # Create the output list.
    my @retVal;
    # Get the type hash.
    my $typeHash = $self->{typeHash};
    # It is now time to create the HTML lines. Start with a heading.
    push @retVal, CGI::div({ class => 'heading' }, CGI::h1($modName));
    push @retVal, CGI::start_div({ id => 'Pod' });
    # Put in the top-of-page anchor.
    push @retVal, CGI::a({ class => 'dummyTopAnchor', name => '___top' });
    # Next comes the table of contents.
    push @retVal, CGI::div({ class => 'indexgroup' }, CGI::ul({ class => 'indexlist indexlist1' },
            CGI::li({ class => 'indexItem indexItem1' }, [ map { CGI::a({ href => "#$_"}, $_) } sort keys %$typeHash ])));
    # Now we loop through the types. The tricky part here is that we need to look inside the type text for
    # callbacks to other types.
    for my $type (sort keys %$typeHash) {
        my $definition = $typeHash->{$type};
        my ($structure, $comment) = @$definition;
        # Generate the section heading.
        push @retVal, CGI::h1(CGI::a({ class => 'u', href => '#___top', name => $type }, $type));
        # Add the comment (if any).
        if ($comment) {
            push @retVal, CGI::p($comment);
        }
        # Process according to the definition type.
        my $html = $self->DisplayDefinition($structure);
        push @retVal, $html;
    }
    # Close the section.
    push @retVal, CGI::end_div(), CGI::br({ class => 'clear' });
    return join("\n", @retVal);
}

=head3 DisplayDefinition

    my $html = $parser->DisplayDefinition($definition);

Compute the HTML for a type definition.

=over 4

=item definition

A type definition, which is either a scalar (indicating another type), or a list consisting of a keyword (C<structure>, C<tuple>, C<mapping>, C<list>)
and a list of elements. Each element is a 3-tuple consisting of (0) a type definition, (1) a comment, and (2) an optional name.

=item RETURN

Returns the HTML to display the type.

=cut

sub DisplayDefinition {
    my ($self, $definition) = @_;
    my $retVal;
    if (! ref $definition) {
        if ($self->{typeHash}{$definition}) {
            $retVal = CGI::a({ href => "#$definition" }, $definition);
        } else {
            $retVal = $definition;
        }
    } else {
        my ($category, @elements) = @$definition;
        if ($category eq 'list') {
            $retVal = CGI::p("list of " . $self->DisplayDefinition($elements[0][0]));
        } elsif ($category eq 'mapping') {
            $retVal = join("\n", CGI::start_table(),
                    CGI::Tr(CGI::th({ colspan => 2 }, '<b><center>MAPPING</center></b>')),
                    CGI::Tr(CGI::td('key'), CGI::td($self->DisplayDefinition($elements[0][0]))),
                    CGI::Tr(CGI::td('value'), CGI::td($self->DisplayDefinition($elements[1][0]))),
                    CGI::end_table());
        } elsif ($category eq 'tuple' || $category eq 'structure') {
            my @table = (CGI::start_table(-title => $category),
                    CGI::Tr((CGI::th({ colspan => 3 }, '<b><center>' . uc($category) . '</center></b>'))));
            for my $element (@elements) {
                my ($type, $comment, $name) = @$element;
                push @table, CGI::Tr(CGI::td($name), CGI::td($self->DisplayDefinition($type)), CGI::td($comment));
            }
            push @table, CGI::end_table();
            $retVal = join("\n", @table);
        }
    }
    return $retVal;
}


=head3 ParseTokens

    ParseTokens(\@tokens)

Consume a single type definition and store it in the type hash.

=over 4

=item tokens

Reference to a list of tokens.

=back

=cut

sub ParseTokens {
    my ($self, $tokens) = @_;
    # Get the next non-comment token.
    my $next = $self->next($tokens);
    if ($next eq 'funcdef') {
        $self->ParseFunction($tokens);
    } elsif ($next eq 'typedef') {
        $self->ParseType($tokens);
    } else {
        die "Invalid token '$next' encountered.";
    }
}

=head3 ParseFunction

    $parser->ParseFunction($tokens);

Consume a function definition. We don't keep these.

=over 4

=item tokens

Reference to a list of tokens.

=back

=cut

sub ParseFunction {
    my ($self, $tokens) = @_;
    # Loop until we reach the end.
    my $done;
    while (! $done) {
        my $token = $self->next($tokens);
        if (! $token || $token eq ';') {
            $done = 1;
        }
    }
    # Clear the comment.
    $self->{comment} = '';
}


=head3 ParseType

    $parser->ParseType($tokens);

Parse and consume a type definition. The type definition will be added to the type hash.

=over 4

=item tokens

Reference to a list of tokens.

=back

=cut

sub ParseType {
    my ($self, $tokens) = @_;
    # Get the definition.
    my ($name, $description, $comment) = $self->ParseDefinition($tokens);
    # Store it in the type hash.
    $self->{typeHash}{$name} = [$description, $comment];
}


=head3 ParseDefinition

    my ($description, $name) = $parser->ParseDefinition($tokens);

Consume and parse a definition. A definition contains a type description followed by a name and then a semi-colon.

=over 4

=item tokens

Reference to a list of tokens.

=item RETURN

Returns a 3-element list consisting of a type description, name, and comment.

=back

=cut

sub ParseDefinition {
    my ($self, $tokens) = @_;
    # Save the comment.
    my $comment = $self->GetComment();
    # Get the type description.
    my $description = $self->ParseDescription($tokens);
    # Get the name.
    my $name = $self->next($tokens);
    # Consume the semi-colon.
    my $delim = $self->next($tokens);
    die "Expected ';' but found '$delim' in definition of $name." if $delim ne ';';
    # Clear any leftover comments.
    $self->{comment} = '';
    # Return the results.
    return ($name, $description, $comment);
}


=head3 ParseDescription

    my $description = $parser->ParseDescription($tokens);

Consume and format a type description. A type description can be a type name, or it can be a complex type-- structure, mapping, list, or tuple. A
structure contains named members each of the form

    description name;

enclosed in braces. A mapping, list, or tuple contains one or more descriptions separated by commas and enclosed in angle brackets.

=over 4

=item tokens

Reference to a list of tokens.

=item RETURN

Returns a list reference describing the type.

=back

=cut

use constant TYPED => {mapping => 1, list => 1, tuple => 1};

sub ParseDescription {
    my ($self, $tokens) = @_;
    # This will be the description to return.
    my $retVal;
    # Get the first token.
    my $category = $self->next($tokens);
    if ($category eq 'structure') {
        $self->ParseDelim($tokens);
        my @description;
        # Here we have a structure with members. Consume the structure.
        while ($tokens->[0] ne '}') {
            my ($name, $element, $comment) = $self->ParseDefinition($tokens);
            push @description, [$element, $comment, $name];
        }
        $retVal = [structure => @description];
        # Consume the trailing delimiter.
        $self->next($tokens);
    } elsif (TYPED->{$category}) {
        $self->ParseDelim($tokens);
        my @description;
        # Here we have a tuple with members. Consume the tuple.
        my $last = '';
        while ($last ne '>') {
            my $element = $self->ParseDescription($tokens);
            my $name = '';
            if ($tokens->[0] ne ',' && $tokens->[0] ne '>') {
                $name = $self->next($tokens);
            }
            # This will be either a comma or end-of-tuple.
            $last = $self->next($tokens);
            push @description, [$element, '', $name];
        }
        $retVal = [$category => @description];
    } else {
        # Here we have a redefined type.
        $retVal = $category;
    }
    # Return the result.
    return $retVal;
}

=head3 ParseDelim

    $parser->ParseDelim($tokens);

Consume a delimiter and update the delimiter stack.

=over 4

=item tokens

Reference to a list of tokens.

=back

=cut

sub ParseDelim {
    my ($self, $tokens) = @_;
    my $delim = $self->next($tokens);
    die "Expected { or <, found closing delimiter or end-of-file." if (! $delim);
}


=head3 GetComment

    my $comment = $parser->GetComment();

Retrieve the last comment and clear it out so that it is not re-used.

=cut

sub GetComment {
    my ($self) = @_;
    my $retVal = $self->{comment};
    $self->{comment} = '';
    return $retVal;
}


=head3 next

    my $token = $parser->next($tokens);

Return the next non-comment token. The last comment found will be stored in this object and the consumed tokens will be pulled off the list.

=over 4

=item tokens

Reference to a list of unconsumed tokens.

=item RETURN

Returns the next non-comment token, or C<undef> if we are at the end of the current section.

=back

=cut

sub next {
    my ($self, $tokens) = @_;
    my $retVal;
    my $done;
    while (@$tokens && ! $done) {
        my $token = shift @$tokens;
        if (index($token, '/*') >= 0) {
            my $comment = $token =~ /\/\*(.+)\*\//s;
            $self->{comment} = $1;
        } else {
            $retVal = $token;
            $done = 1;
        }
    }
    return $retVal;
}


=head3 ShowType

    my $html = ShowType($type, $typeHash);

Display a type name in HTML. If the name is found in the type hash, we
wrap it in a link.

=over 4

=item type

Type name to display.

=item typeHash

Reference to a hash mapping each other type to its definition. This is used to determine
if we can link a type name to its section in the HTML document.

=item RETURN

Returns the type name in HTML.

=back

=cut

sub ShowType {
    # Get the parameters.
    my ($type, $typeHash) = @_;
    # Declare the return variable.
    my $retVal = $type;
    # Analyze the type.
    if ($type =~ /list\s*<(.+)>/) {
        # Here we have a list of another type.
        $retVal = "list of " . ShowType($1, $typeHash);
    } elsif ($type =~ /tuple\s*<(.+)>/) {
        # Here we have a tuple of other types.
        my @parts = split /\s*,\s*/, $1;
        $retVal = 'tuple of (' .  join(", ", map { ShowType($_, $typeHash) } @parts) . ')';
    } elsif (exists $typeHash->{$type}) {
        # Here we have a linkable type.
        $retVal = CGI::a({ href => "#$type" }, $retVal);
    } elsif ($type =~ /mapping\s*<(.+),(.+)>/) {
        # Here we have a hash map.
        my ($source, $target) = ($1, $2);
        $retVal = 'mapping of ' . ShowType($source, $typeHash) . ' to ' . ShowType($target, $typeHash);
    } elsif ($type =~ /structure\s+{(.+)}/) {
        # Here we have a sub-structure.
        my @lines = StructureBody($1, $typeHash, ';');
        $retVal = join("\n", @lines);
    }
    # Return the result.
    return $retVal;
}


=head3 FormatComment

    my $html = FormatComment($text);

Format a comment for HTML. This involves looking for special indicators
that the text is not simply raw word sequences.

=over 4

=item text

Comment text to convert.

=item RETURN

A version of the comment text suitable for the interior of an HTML paragrapph.

=back

=cut

sub FormatComment {
    # Get the parameters.
    my ($text) = @_;
    # Declare the return variable.
    my $retVal = $text;
    # Remove excess spaces.
    $retVal =~ s/^\*\s+//;
    $retVal =~ s/\n\s+\*/ /gs;
    # Return the result.
    return CGI::i($retVal);
}

1;
