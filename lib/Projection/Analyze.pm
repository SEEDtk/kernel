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


package Projection::Analyze;

    use strict;
    use warnings;
    use CGI;
    use SEEDClient;
    use SeedUtils;

## TODO: output the json of this object. Create a CGI to verify against core seed and only produce stuff still relevant (genome not in subsystem, peg has same function)

=head1 Methods for Analyzing Projection Output

This object analyzes the output of a projection run (see L<solid_projection_engine.pl>). It will read the
C<solid.proj.tbl> and C<core.proj.tbl> files and produce two hashes-- one that mimics a spreadsheet of the
projected genomes and one that lists the genomes that appear to be inactive in the subsystem but have suspicious
features possessing the subsystem roles.

=head2 Object Description

=head3 Features

In the hashes, each feature is described by a sub-hash consisting of the following fields

=over 4

=item id

The feature ID.

=item fun

The function ID.

=item status

C<good> if the feature is good; otherwise a description of why the feature is suspicious

=item role

The ID of the role performed in the subsystem.

=item check

The feature's checksum, indicating its status when loaded.

=back

=head3 Spreadsheet

In the spreadsheet hash, each genome ID maps to a hash containing the following fields.

=over 4

=item variant

The variant code.

=item name

The genome name.

=item model

The model genome used to predict the variant.

=item status

C<new> if this is a new projection, C<changed> if it is a projection that matches an existing database entry but has
different features, and C<old> if it is a projection that exactly matches an existing database entry.

=item row

Reference to a list containing the spreadsheet cells, in order. Each cell in turn contains a list of feature hashes.

=back

=head3 Object Fields

The fields of this object are as follows.

=over 4

=item subsystem

The subsystem ID.

=item name

The subsystem name.

=item roles

Reference to a list of roles for the subsystem, in order.

=item roleMap

Reference to a hash mapping each role ID to a 2-tuple consisting of (0) the role's
ordinal position in the subsystem, and (1) the role description.

=item db

The L<Shrub> object for accessing the database.

=item sheet

Reference to a hash containing the subsystem spreadsheet with the proposed changes.

=item errors

Reference to a hash keyed on genome ID that maps each genome to a list of suspicious features.

=item genomes

Reference to a hash of genome IDs to names for all the acceptable genomes.

=item seedRoot

URL of the seedviewer page for HTML links.

=item core

TRUE if only core genoems are being processed, else FALSE.

=item pegMap

Reference to a hash mapping each genome ID to a sub-hash that maps its peg IDs to their checksums.

=back

=head2 Special Methods

=head3 new

    my $proj = Projection::Analyze->new($shrub, $subID, %options);

Create a new, blank subsystem projection analysis object.

=over 4

=item shrub

The L<Shrub> object for accessing the database.

=item subID

The ID of the relevant subsystem.

=item options

A hash of options which may include zero or more of the following keys.

=over 8

=item core

If TRUE, then only core genome projections will be included in the output. The default is FALSE.

=back

=back

=cut

sub new {
    # Get the parameters.
    my ($class, $shrub, $subID, %options) = @_;
    # Initialize the object.
    my $retVal = {
        subsystem => $subID,
        db => $shrub,
        sheet => {},
        errors => {},
        pegMap => {},
        core => ($options{core} // 0)
    };
    # Get the subsystem name.
    my ($name) = $shrub->GetFlat('Subsystem', 'Subsystem(id) = ?', [$subID], 'name');
    # Get the roles.
    my @roles = $shrub->GetAll('Subsystem2Role Role', 'Subsystem2Role(from-link) = ? ORDER BY Subsystem2Role(ordinal)',
            [$subID], 'to-link ordinal Role(description)');
    # Store them in the object.
    $retVal->{name} = $name;
    $retVal->{roles} = [map { $_->[0] } @roles];
    # Compute the role map.
    my %roleMap = map { $_->[0] => [$_->[1], $_->[2]] } @roles;
    $retVal->{roleMap} = \%roleMap;
    # Create the genome map. Only the genomes we are going to allow in the output will be included.
    my ($filter, $parms);
    if ($options{core}) {
        $filter = 'Genome(core) = ?';
        $parms = [1];
    } else {
        $filter = '';
        $parms = [];
    }
    my %genomes = map { $_->[0] => $_->[1] } $shrub->GetAll('Genome', $filter, $parms, 'id name');
    $retVal->{genomes} = \%genomes;
    # Bless and return the object.
    bless $retVal, $class;
    return $retVal;
}

=head3 new_from_file

    my $proj = Projection::Analyze->new_from_file($shrub, $fileName);

Read this object from a file produced by L</JsonOut>.

=over 4

=item shrub

L<Shrub> object for accessing the database.

=item fileName

Name of a file containing this object in JSON format.

=back

=cut

sub new_from_file {
    # Get the parameters.
    my ($class, $shrub, $fileName) = @_;
    # Read the object.
    my $retVal = SeedUtils::read_encoded_object($fileName);
    # Add the database.
    $retVal->{db} = $shrub;
    # Bless and return the object.
    bless $retVal, $class;
    return $retVal;
}


=head2 Public Manipulation Methods

=head3 JsonOut

    $proj->JsonOut($fileName);

Write this object out to a file in JSON format.

=over 4

=item fileName

Name of the file to receive the object in JSON format.

=back

=cut

sub JsonOut {
    # Get the parameters.
    my ($self, $fileName) = @_;
    # We need to create an unblessed object to hold all our data. We copy all the fields except "db".
    my %output;
    for my $key (keys %$self) {
        if ($key ne 'db') {
            $output{$key} = $self->{$key};
        }
    }
    # Write the object.
    SeedUtils::write_encoded_object(\%output, $fileName);
}

=head3 AddFile

    $proj->AddFile($fileName);

Merge the data from the specified file into this object. The file must be one of the output files
(C<solid.proj.tbl> or C<core.proj.tbl>) from L<solid_projection_engine.pl>. Any projections already
in the object will be overwritten if there is a conflict. (That is, if data in the file is for a genome
that is already in the object, the old genome's data will be erased.)

=over 4

=item fileName

Name of the file to process.

=back

=cut

sub AddFile {
    # Get the parameters.
    my ($self, $fileName) = @_;
    # Open the input file.
    open(my $ih, "<", $fileName) || die "Could not open projection file: $!";
    # Loop through the input file, reading prediction sections.
    while (! eof $ih) {
        my ($genome, $vc, $model, $cells, $errors) = $self->_ReadPrediction($ih);
        # This will hold the genome's peg list.
        my %pegList;
        # Only proceed if we care about this genome.
        my $gName = $self->{genomes}{$genome};
        if ($gName) {
            if ($vc) {
                # Here we have a prediction of an active variant. Check the database.
                my ($vcOld) = $self->{db}->GetFlat('Genome2Row SubsystemRow Row2Subsystem',
                        'Genome2Row(from-link) = ? AND Row2Subsystem(to-link) = ? ORDER BY SubsystemRow(privilege) DESC',
                        [$genome, $self->{subsystem}], 'SubsystemRow(variant-code)');
                my $status = 'old';
                if (! $vcOld) {
                    $status = 'new';
                } elsif ($vcOld ne $vc) {
                    $status = 'changed';
                }
                $self->{sheet}{$genome} = { variant => $vc, name => $gName, model => $model, status => $status, row => $cells };
                # Get the list of peg IDs for the genome.
                for my $cell (@$cells) {
                    for my $peg (@$cell) {
                        $pegList{$peg->{id}} = $peg->{check};
                    }

                }
            }
            # If we have errors, store them.
            if (@$errors) {
                $self->{errors}{$genome} = $errors;
                # Update the list of peg IDs.
                for my $error (@$errors) {
                    $pegList{$error->{id}} = $error->{check};
                }
            }
            # Store the peg ID list.
            $self->{pegMap}{$genome} = [sort keys %pegList];
        }
    }
}

=head3 ToHtml

    $proj->ToHtml($fileName, $seedRoot);

Write the analysis to the specified output file in HTML table form, with links to the specified seed.

=over 4

=item fileName

Name of the HTML output file.

=item seedRoot

Root URL for seed viewer pages in the target SEED (used in hyperlinks).

=back

=cut

use constant GENOME_COLOR => { 'new' => '#88FF88', 'changed' => '#FFFF00' };
use constant FEATURE_ERROR => '#FFFF00';
use constant TABLE_DOUBLED => '#FFCCCC';

sub ToHtml {
    # Get the parameters.
    my ($self, $fileName, $seedRoot) = @_;
    # Save the seed root URL.
    $self->{seedRoot} = $seedRoot;
    # If this is core-only mode, get the seed client object.
    my $seedClient;
    if ($self->{core}) {
        $seedClient = SEEDClient->new('http://core.theseed.org/FIG/seed_svc');
    }
    # Open the output file.
    open(my $oh, ">", $fileName) || die "Could not open HTML output file: $!";
    # Start the page.
    my $title = "$self->{name} Projection Analysis";
    _Print($oh, CGI::start_html(-title => $title, -target => '_blank'),
        '<style>',
        'table, th, td, tr { border-style: double; border-collapse: collapse; }',
        'th { text-align: left; }',
        '</style>',
        CGI::h1($title));
    # Start the projections table. Note the heading columns. We have the genome,
    # the variant code, the model genome ID, and a list of role IDs.
    _Print($oh, CGI::h2("Projections"), CGI::start_table(),
            CGI::Tr(CGI::th("Genome"), CGI::th("VC"), CGI::th("Model"),
            map { CGI::th($_) } @{$self->{roles}}));
    # Loop through the projections.
    my $sheet = $self->{sheet};
    for my $genome (sort keys %$sheet) {
        # Get this projection.
        my $projection = $sheet->{$genome};
        # If this is core mode, verify the features.
        my %badPegs;
        if ($self->{core}) {
            my $pegHash = $self->{pegMap}{$genome};
            my @pegList = keys %$pegHash;
            my $funMap = $seedClient->get_function(\@pegList);
            for my $peg (@pegList) {
                my $fun = $funMap->{$peg};
                if (! $fun) {
                    $badPegs{$peg} = 1;
                } else {
                    my $md5 = Shrub::Checksum($fun);
                    if ($md5 ne $pegHash->{$peg}) {
                        $badPegs{$peg} = 1;
                    }
                }
            }
        }
        # Compute the color based on the status. Undefined is a legal result, and
        # indicates no special coloring.
        my $color = GENOME_COLOR->{$projection->{status}};
        # Format the genome and name.
        my $genomeHead = $self->_Genome($genome, $projection->{name}, $color);
        # Loop through the spreadsheet cells, creating table cells.
        my @cellHtmls;
        for my $cell (@{$projection->{row}}) {
            # If the cell is empty, use a non-breaking space.
            if (! @$cell) {
                push @cellHtmls, CGI::td('&nbsp;');
            } else {
                # Color the cell if there are multiple pegs in this cell.
                my %attr;
                if (@$cell > 1) {
                    $attr{style} = 'background-color: ' . TABLE_DOUBLED;
                }
                # Generate the peg data.
                my @pegHtmls;
                for my $feature (@$cell) {
                    my $id = $feature->{id};
                    my ($title, $color);
                    if ($feature->{status} ne 'good') {
                        $title = $feature->{status};
                        $color = FEATURE_ERROR;
                    }
                    push @pegHtmls, $self->_Feature($id, $title, $color);
                }
                push @cellHtmls, CGI::td(\%attr, join(", ", @pegHtmls));
            }
        }
        # Generate the table row.
        _Print($oh, CGI::Tr(CGI::th($genomeHead), CGI::td($projection->{variant}),
                CGI::td($self->_Genome($projection->{model}), @cellHtmls)));
    }
    # Close the projections table.
    _Print($oh, CGI::end_table());
    # Start the error table.
    _Print($oh, CGI::h2('Suspicious Features'), CGI::start_table(),
            CGI::Tr(CGI::th('Genome'), CGI::th('Feature'), CGI::th('Error')));
    # Loop through the genomes.
    my $errors = $self->{errors};
    for my $genome (sort keys %$errors) {
        # Get this genome's error list.
        my $errorList = $errors->{$genome};
        # Format the genome HTML.
        my $genomeHtml = $self->_Genome($genome, $self->{genomes}{$genome});
        # Loop through the error features.
        for my $feature (@$errorList) {
            # Display this feature.
            _Print($oh, CGI::Tr(CGI::th($genomeHtml), CGI::td($self->_Feature($feature->{id}),
                    CGI::td($feature->{status}))));
        }
    }
    # Close the error table.
    _Print($oh, CGI::end_table());
    # Close the page.
    _Print($oh, CGI::end_html());
}


=head2 Internal Utility Methods

=head3 _Genome

    my $html = $proj->_Genome($genomeID, $name, $color);

Format the HTML to display a genome link. This includes linking the genome ID to the genome page.

=over 4

=item genomeID

ID of the genome to which we are linking.

=item name (optional)

If specified, the name of the genome, to be included in the link text.

=item color (optional)

The background color to give to the link.

=item RETURN

Returns the HTML for a genome link, titled with the genome name, displaying the genome ID and
optionally the genome name, linking to the genome page in the stored seed root.

=back

=cut

sub _Genome {
    # Get the parameters.
    my ($self, $genomeID, $name, $color) = @_;
    # Format the genome text.
    my $text = $genomeID;
    if ($name) {
        $text .= " $name";
    }
    # Compute the title.
    my $title = $self->{genomes}{$genomeID};
    # Compute the link attributes.
    my %attrs = (href => "$self->{seedRoot}?page=Organism;organism=$genomeID", title => $title);
    # Compute the style.
    if (defined $color) {
        $attrs{style} = "background-color: $color";
    }
    # Create the link HTML.
    my $retVal = CGI::a(\%attrs, $text);
    # Return the result.
    return $retVal;
}

=head3 _Feature

    my $html = $proj->_Feature($fid, $title, $color);

Format the HTML to display a feature. This includes linking the feature ID to the SEED feature page.

=over 4

=item fid

ID of the feature to display.

=item title

Title to display when hovering over the feature ID.

=item color (optional)

Background color to use for the link.

=item RETURN

Returns the HTML for a link to the feature page, titled with the specified string, and optionally
colored with the specified hue.

=back

=cut

sub _Feature {
    # Get the parameters.
    my ($self, $fid, $title, $color) = @_;
    # Compute the link attributes.
    my %attr = (href => "$self->{seedRoot}?page=Annotation;feature=$fid", title => $title);
    if (defined $color) {
        $attr{style} = "background-color: $color";
    }
    # Adjust the feature ID by removing the genome ID and defaulting to pegs.
    my $adjusted = $fid;
    if ($fid =~ /\.peg\.(\d+)$/) {
        $adjusted = $1;
    } elsif ($fid =~ /\.(\w+\.\d+)$/) {
        $adjusted = $1;
    }
    # Create the link HTML.
    my $retVal = CGI::a(\%attr, $adjusted);
    # Return the result.
    return $retVal;
}

=head3 _ReadPrediction

    my ($genome, $vc, $model, $pegs, $errors) = $proj->_ReadPrediction($ih);

Read and analyze a subsystem prediction from the specified input stream.
The subsystem prediction consists of a header record that contains the
genome ID, variant code, and model genome ID, zero or more peg
predictions, and zero or more supicious-peg (error) statements. The exact
format is described in L<solid_projection_engine.pl/Output>.

=over 4

=item ih

Open handle to an input stream positioned at the beginning of a prediction.

=item RETURN

Returns a 5-element list consisting of the following elements.

=over 8

=item genome

The ID of the genome containing the predicted subsystem instance.

=item vc

The variant code, or an empty string if the subsystem is not active in this genome.

=item model

The ID of the genome used as a model for predicting the variant.

=item pegs

Reference to a list containing one entry per subsystem role. Each list entry will
contain a sub-list of all the L</Features> performing the role.

=item errors

Reference to a list of L</Features> that have suspicious characteristics.

=back

=back

=cut

use constant BAD_VARIANTS => { 'dirty' => 1, 'not-active' => 1 };

sub _ReadPrediction {
    # Get the parameters.
    my ($self, $ih) = @_;
    # Declare the error return list.
    my @errors;
    # Initialize the spreadsheet peg list.
    my @cells = map { [] } @{$self->{roles}};
    # Get the header record.
    my ($subID, $genome, $vc, $model) = _ReadRecord($ih);
    # We will put the pegs we find in here. Each peg will map to a feature object.
    my %pegH;
    # If this is a bad variant, skip it.
    if (BAD_VARIANTS->{$vc}) {
        # Loop until we hit the separator. This spins until we read a line that begins with
        # four hyphens.
        while (! eof $ih && <$ih> !~ /^----/) {}
        # Denote this is a bad variant.
        $vc = '';
    } else {
        # Loop through the pegs.
        my $done;
        while (! eof $ih && ! $done) {
            my ($flag, $peg, $fun, $role, $checksum) = _ReadRecord($ih);
            # Check for the end marker.
            if ($flag eq '----') {
                $done = 1;
            } else {
                # Store this peg.
                $pegH{$peg} = { id => $peg, role => $role, fun => $fun, status => 'good', check => $checksum };
            }
        }
    }
    # Loop through the errors.
    my $done = 0;
    while (! eof $ih && ! $done) {
        my ($peg, $error, $desc, $fun, $role, $checksum) = _ReadRecord($ih);
        # Check for the end marker.
        if ($peg eq '//') {
            $done = 1;
        } else {
            # Store this error.
            if ($desc) {
                $error .= " ($desc)";
            }
            my $pegData = { id => $peg, role => $role, fun => $fun, status => $error, check => $checksum };
            push @errors, $pegData;
            # Add it to the subsystem if the subsystem is being kept.
            if ($vc) {
                $pegH{$peg} = $pegData;
            }
        }
    }
    # Check for wrong pegs.
    # If we have spreadsheet pegs, format the spreadsheet.
    my $roleMap = $self->{roleMap};
    for my $peg (keys %pegH) {
        my $pegData = $pegH{$peg};
        my $role = $pegData->{role};
        push @{$cells[$roleMap->{$role}[0]]}, $pegData;
    }
    # Return the results.
    return ($genome, $vc, $model, \@cells, \@errors);
}

=head3 _ReadRecord

    my @fields = _ReadRecord($ih);

Read a record from the input and separate it into tab-delimited fields.

=over 4

=item ih

Open input stream from which to read a record.

=item RETURN

Returns a list of the tab-delimited fields in the record.

=back

=cut

sub _ReadRecord {
    # Get the parameters.
    my ($ih) = @_;
    # Read the line.
    my $line = <$ih>;
    chomp $line;
    # Declare the return variable.
    my @retVal = split /\t/, $line;
    # Return the result.
    return @retVal;
}

=head3 _Print

    _Print($oh, @lines);

Print the specified lines, with separating new-lines.

=over 4

=item oh

Open output handle for the output file.

=item lines

List of lines to print.

=back

=cut

sub _Print {
    my ($oh, @lines) = @_;
    for my $line (@lines) {
        print $oh "$line\n";
    }
}

1;