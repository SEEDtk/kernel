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


package gjo::BlastParse;

    use strict;
    use warnings;
    use Carp;

=head1 Blast Parser

This object parses BLAST+ output. It is an object-oriented replacement for the old C<gjoparseblast>
module. This module used a clever trick to handle input state information that should only have
failed on extremely rare occasions; however, our clever and talented staff found a way to make it
fail every time. This replacement module uses object-oriented techniques to handle the state
information, which is more reliable. It does not yet have the full powers of its predecessor.

The basic usage strategy is

    my $parser = gjo::BlastParse->new($input, %options);
    while (my $hsp = $parser->next_hsp) {
        ... process ...
    }

An hsp is a tuple of the following form

     [ qid qdef qlen sid sdef slen scr e_val p_n p_val n_mat n_id n_pos n_gap dir q1 q2 qseq s1 s2 sseq ]
        0   1    2    3   4    5    6    7    8    9    10    11   12    13   14  15 16  17  18 19  20

Alternatively, you can use

    my $parser = gjo::BlastParse->new($input, %options);
    while (my $rec = $parser->next_record) {
        ... process ...
    }

where a record is one of the following

     [ 'Query='  query_id  query_def  query_len ]
          0         1          2          3

     [   '>'     sbjct_id  sbjct_def  sbjct_len ]
          0         1          2          3

     [ 'HSP' scr exp p_n p_val n_mat n_id n_sim n_gap dir q1 q2 qseq s1 s2 sseq ]
         0    1   2   3    4     5    6     7     8    9  10 11  12  13 14  15

You can also ask for all the records in the form of a list.

    my @hsps = gjo::BlastParse::all_hsps($input, %options);

or

    my @records = gjo::BlastParse::all_records($input, %options);


This object contains the following fields.

=over 4

=item input

Input source for this parsing process. Either an array reference or a file input handle.

=item inpos

If the input is an array reference, the index into the array of the next record; otherwise
undefined.

=item query

Current query record in the blast output.

=item subject

Current subject record in the blast output.

=item options

Reference to the options hash (see L</new>).

=item line

Input line to process.

=back

=head2 Special Methods

=head3 new

    my $varname = gjo::BlastParse->new($input, %options);

Create a new Blast Parser for a particular input stream.

=over 4

=item input

Input stream containing the blast output. Either an open input file handle, a string broken by new-lines, or
a reference to an array. If C<undef>, the standard input will be used.

=item options

A hash of options containing zero or more of the following keys.

=over 8

=item self

If TRUE, then matches between a sequence and itself will be included in the output. The default
is FALSE (self-matches will be discarded).

=back

=back

=cut

sub new {
    # Get the parameters.
    my ($class, $input, %options) = @_;
    # Determine the input.
    my ($inputField, $pos);
    if (! defined $input) {
        # No input specified, use STDIN.
        $inputField = \*STDIN;
    } elsif (! ref $input) {
        # A scalar. Convert it into an array.
        $inputField = [ split "\n", $input ];
        $pos = 0;
    } elsif (ref $input eq 'GLOB') {
        # Open file handle.
        $inputField = $input;
    } elsif (ref $input eq 'ARRAY') {
        # Array of data lines.
        $inputField = $input;
        $pos = 0;
    } else {
        confess("Invalid input type " . ref($input) . " specified for BlastParse.");
    }
    # Create the object. Note that the query and subject start out blank.
    my $retVal = {
        input => $inputField,
        inpos => $pos,
        query => ['Query=', '', '', 0],
        subject => ['>', '', '', 0],
        options => { %options }
    };
    # Get the first input line.
    _nextline($retVal);
    # Bless and return it.
    bless $retVal, $class;
    return $retVal;
}


=head2 List Methods

=head3 all_hsps

    my @hsps = gjo::BlastParse::all_hsps($input, %options);

Return the blast output in hsp format.

=over 4

=item input

Input stream containing the blast output. Either an open input file handle, a string broken by new-lines, or
a reference to an array. If C<undef>, the standard input will be used.

=item options

A hash of options containing zero or more of the following keys.

=over 8

=item self

If TRUE, then matches between a sequence and itself will be included in the output. The default
is FALSE (self-matches will be discarded).

=back

=item RETURN

Returns a list of 21-tuples containing the blast output in HSP format. Each tuple contains (0) the
query ID, (1) the query description, (2) the query sequence length, (3) the subject ID, (4) the
subject description, (5) the subject length, (6) the match score, (7) the e-value, (8) the p-number,
(9) the Poisson score, (10) the match length, (11) the identity count, (12) the positive count,
(13) the gap count, (14) the frame direction or shift, (16) the query end position, (17) the query
alignment sequence, (18) the subject start position, (19) the subject end position, and (20) the subject sequence.

=back

=cut

sub all_hsps {
    my ($input, %options) = @_;
    my $helper = gjo::BlastParse->new($input, %options);
    my @retVal;
    while (my $hsp = $helper->next_hsp) {
        push @retVal, $hsp;
    }
    return @retVal;
}

=head3 all_records

    my @records = gjo::BlastParse::all_records($input, %options);

Return the blast output as structured records. Each record is either q query, a subject, or an hsp.

=over 4

=item input

Input stream containing the blast output. Either an open input file handle, a string broken by new-lines, or
a reference to an array. If C<undef>, the standard input will be used.

=item options

A hash of options containing zero or more of the following keys.

=over 8

=item self

If TRUE, then matches between a sequence and itself will be included in the output. The default
is FALSE (self-matches will be discarded).

=back

=item RETURN

Returns a list of queries, subjects, and hsps.

=over 8

=item 1

Each query is a 4-tuple consisting of (0) the string C<Query=>, (1) the query ID, (2) the
query description, and (3) the query sequence length.

=item 2

Each subject is a 4-tuple consisting of (0) the string C<< > >>, (1) the subject ID, (2) the
subject description, and (3) the subject sequence length.

=item 3

Each HSP is a 16-tuple containing (0) the string C<HSP>, (1) the match score, (2) the e-value,
(3) the p-number, (4) the Poisson value, (5) the match length, (6) the identity count,
(7) the positive count, (8) the gap count, (9) the frame direction or shift, (10) the query
start position, (11) the query end position, (12) the query alignment sequence,
(13) the subject start position, (14) the subject end position, and (15) the subject sequence.

=back

=back

=cut

sub all_records {
    my ($input, %options) = @_;
    my $helper = gjo::BlastParse->new($input, %options);
    my @retVal;
    while (my $rec = $helper->next_record) {
        push @retVal, $rec;
    }
    return @retVal;
}

=head2 Iterative Methods

=head3 next_record

    my $record = $parser->next_record

Return the next data record for the output. The data record will either be a query specification, a subject
specification, or an HSP (see L</all_records>).

=cut

sub next_record {
    # Get the parameters.
    my ($self) = @_;
    # This will contain the return value.
    my $retVal;
    # Get the current line.
    my $line = $self->{line};
    # Loop until we run out of input or find a heading line.
    while (! $retVal && defined $line) {
        if ($line =~ /^(Query=|>)\s*(\S.*?)\s*$/) {
            $retVal = $self->_processSequenceRecord($1, $2);
            # We must save sequence records.
            if ($retVal->[0] eq 'Query=') {
                $self->{query} = $retVal;
            } else {
                $self->{subject} = $retVal;
            }
        } elsif ($line =~ /^ Score =\s+([\d.e+-]+.*?)\s*$/) {
            # Here we have an HSP. Do we want to keep it?
            if ($self->{options}{self} || $self->{query}[1] ne $self->{subject}[1]) {
                # Yes. Process it.
                $retVal = $self->_processHSP($1);
            } else {
                # No. Here we have a self-match and we are skipping those.
                while (defined $line) {
                    $line = $self->_next_hsp_line;
                }
                # Try again.
                $line = $self->{line};
                $retVal = undef;
            }
        } else {
            # Nothing special. Keep looking.
            $line = $self->_nextline;
        }
    }
    # Return the record found.
    return $retVal;
}

=head3 next_hsp

    my $record = $parser->next_hsp;

Return an HSP record for the output. The HSP will contain a query specification, a subject
specification, and the HSP data (see L</all_hsps>).

=cut

sub next_hsp {
    # Get the parameters.
    my ($self) = @_;
    # Loop through records until we find an HSP.
    my $retVal = $self->next_record;
    while ($retVal && $retVal->[0] ne 'HSP') {
        $retVal = $self->next_record;
    }
    if ($retVal) {
        # We have an HSP. Pop off the type indicator.
        shift @$retVal;
        # Add the query and subject.
        $retVal = [@{$self->{query}}[1,2,3], @{$self->{subject}}[1,2,3], @$retVal];
    }
    # Return the HSP record.
    return $retVal;
}

## TODO BlastParse all_queries [Q, [S, [H, H ...]], [S, [H, H, ...]], ...], ...
## TODO BlastParse all_subjects [Q, S, [H, H, ...]], ...

=head2 Internal Utility Methods

=head3 _nextline

    my $line = $parser->_nextline()

Return the next line of input.

=cut

sub _nextline {
    # Get the parameters.
    my ($self) = @_;
    # The return value goes in here.
    my $retVal;
    # Check our position in the input.
    my $inpos = $self->{inpos};
    if (defined $inpos) {
        # Here we are cycling through an array.
        $retVal = $self->{input}[$inpos];
        $self->{inpos} = $inpos + 1;
    } else {
        # Here we are reading a file handle.
        my $ih = $self->{input};
        $retVal = <$ih>;
    }
    # Save the line read and return it.
    $self->{line} = $retVal;
    return $retVal;
}

=head3 _processSequenceRecord

    my $record = $self->_processSequenceRecord($type, $data);

Return a query or subject record that begins on the current line of the input
stream. We must extract the ID, description, and length. When done,
the input stream will be positioned on the length line.

=over 4

=item type

The type of sequence record-- C<Query=> for a query, C<< > >> for a subject.

=item data

Text on the query header line after the query marker. It will already be
space-trimmed on both ends.

=item RETURN

Returns a 4-tuple consisting of (0) the string C<Query=>, (1) the query ID, (2) the
query description, and (3) the query sequence length.

=back

=cut

sub _processSequenceRecord {
    # Get the parameters.
    my ($self, $type, $data) = @_;
    # Parse the data. If there is no description, set it to the null string.
    my ($id, $def) = split " ", $data, 2;
    $def //= '';
    # Declare the return variable. We start with type and ID.
    my @retVal = ($type, $id);
    # Look for the length.
    my $length;
    my $line = $self->_nextline;
    while (defined $line && ! defined $length) {
        if ($line =~ /^(?:\s+Length = |Length=)([1-9][\d,]*)$/) {
            # Here we've found the length. Note we have to get rid of the commas.
            $length = $1;
            $length =~ s/,//g;
            # Add the description and length to the return list.
            push @retVal, $def, $length;
        } elsif ($line =~ /^(?:Query=|>|\s+Score\s+=)/) {
            # Here we couldn't find a length.
            confess("No length found for sequence $id in BLAST output.");
        } else {
            # Here we have a continuation of the query description. Trim spaces.
            $line =~ s/^\s+//;
            # If the line is empty, ignore it. If the current description ends in
            # a hyphen, concatenate the new description text, otherwise join with a space.
            if ($line) {
                $def .= ((substr($line, -1, 1) eq '-') ? '' : ' ') . $line;
            }
            # Get the next line.
            $line = $self->_nextline;
        }
    }
    # Return the result.
    return \@retVal;
}


=head3 _processHSP

    my $hspRecord = $parser->_processHSP($data);

Parse BLAST output to produce an HSP record. A great deal of data must be
parsed from the input stream, including the sequence information, which
is split
across several lines in the form of three lines of data per section.

=over 4

=item data

Data from the first line of the HSP section, beginning with the first digit
after the C<Score:> marker.

=item RETURN

Returns a 16-tuple containing (0) the string C<HSP>, (1) the match score, (2) the e-value,
(3) the p-number, (4) the Poisson value, (5) the match length, (6) the identity count,
(7) the positive count, (8) the gap count, (9) the frame direction or shift, (10) the query start
position, (11) the query end position, (12) the query alignment sequence,
(13) the subject start position, (14) the subject end position, and (15) the subject sequence.

=back

=cut

sub _processHSP {
    # Get the parameters.
    my ($self, $data) = @_;
    # Declare the return variable.
    my @retVal = qw(HSP);
    # The various numbers go in here.
    my ($eval, $pn, $pval, $mlen, $ident, $positive, $gap, $frame);
    # Parse the incoming data line.
    my ($score, $remainder) = split " ", $data, 2;
    # Parse the e-value.
    unless ($remainder =~ /Expect(?:\((\d+)\))? =\s+([^\s,]+)/) {
        confess("Invalid Expect syntax--> $data");
    } else {
        $eval = $2;
        if ($1) {
            $pn = $1;
        } elsif ($remainder =~ /P\((\d+)\)/) {
            $pn = $1;
        } else {
            $pn = 1;
        }
    }
    if ($remainder =~ /P(?:\(\d+\))? += +(\S+)$/) {
        $pval = $1;
    } else {
        $pval = $eval;
    }
    # Now we search for the Identities line.
    my $line = $self->_next_hsp_line;
    if (! defined $line) {
        confess("End-of-record searching for identities line--> $data");
    } elsif ($line =~ /^ Identities =\s+(\d+)\/(\d+)/) {
        ($ident, $mlen) = ($1, $2);
    } else {
        confess("Invalid identities line--> $line");
    }
    # Here we extract the other stuff from the Identities line.
    $positive = ($line =~ /Positives =\s+(\d+)/ ? $1 : $ident);
    $gap = ($line =~ /Gaps =\s+(\d+)/ ? $1 : 0);
    # Check for a frame.
    $line = $self->_next_hsp_line;
    if (defined $line && $line =~ /Frame =\s+(\S)/) {
        $frame = $1;
        $line = $self->_next_hsp_line;
    }
    # Finally, we search for the query and subject alignment. This is the
    # hard part. We need to remember the start and stop points for the
    # two sequences and accumulate the sequence strings.
    my %q = (i1 => undef, i2 => 0, seqs => []);
    my %s = (i1 => undef, i2 => 0, seqs => []);
    my %aligns = (Query => \%q, Sbjct => \%s);
    # Loop through the data lines until we hit end-of-record.
    while (defined $line) {
        # If we have an alignment line, process it.
        if ($line =~ /^(Query|Sbjct)\s+(\d+)\s+(\S+)\s+(\d+)/) {
            # Get the hash for the sequence in question.
            my $alignH = $aligns{$1};
            # Extract the start offset, sequence letters, and end offset.
            my ($start, $seq, $end) = ($2, $3, $4);
            # Store the start offset if it is our first time.
            $alignH->{i1} //= $start;
            # Append the sequence.
            push @{$alignH->{seqs}}, $seq;
            # Store the end offset. The last end offset wins.
            $alignH->{i2} = $end;
        }
        # Get the next line.
        $line = $self->_next_hsp_line;
    }
    # Form the result and return it. Start with the simple stuff.
    push @retVal, $score, $eval, $pn, $pval, $mlen, $ident, $positive, $gap;
    # If there is no frame, we must compute it.
    if (! defined $frame) {
        $frame = ((seqdir($q{i1}, $q{i2}) == seqdir($s{i1}, $s{i2})) ? 1 : -1);
    }
    push @retVal, $frame;
    # Add the sequence data.
    push @retVal, $q{i1}, $q{i2}, join("", @{$q{seqs}}), $s{i1}, $s{i2}, join("", @{$s{seqs}});
    # Return the result.
    return \@retVal;
}

=head3 _next_hsp_line;

    my $line = $parser->_next_hsp_line;

Return the next line in the input stream, or C<undef> if we hit the end
of an HSP record. This requires we check for a lot of possible weirdness
in the input line. The line will still be stored internally for the next
operation, but it will not be returned to the caller.

=cut

sub _next_hsp_line {
    # Get the parameters.
    my ($self) = @_;
    # Declare the return variable.
    my $retVal;
    # Get the next line.
    my $line = $self->_nextline;
    # We only return the line if it is defined and DOES NOT match one of the
    # following patterns.
    unless (! defined $line ||
            $line =~ /^(?:Query=|>|Parameters:|Lambda\s+K\s+H\s+a\s+alpha)/ ||
            $line =~ /^\s+(?:Score = |Plus Strand HSPs:|Minus Strand HSPs:|Database:)/) {
        $retVal = $line;
    }
    # Return the result.
    return $retVal;
}

=head3 seqdir

    my $dir = gjo::BlastParse::seqdir($i1, $i2);

Return 1 if the first number is at or to the left of the second, else -1.

=cut

sub seqdir {
    return ( $_[0] <= $_[1] ) ? 1 : -1;
}

1;