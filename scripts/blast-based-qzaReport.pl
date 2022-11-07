#!/usr/bin/env perl

=head1 Driver for Gary Olsen's BLASTN-based Amplicon-to-RepGen analysis toolkit

    blast-based-qzaReport [options] --qza ampliconQzaFile --ssu SsuHypervarRegionsFasta --reportname=qzaReport.tab
or
    blast-based-qzaReport [options] --qza ampliconQzaFile --ssu SsuHypervarRegionsFasta > qzaReport.tab

This script generates a qzaReport-formatted report of read-counts per reference-genome from a QZA-archive of reads.
Output will be written to STDOUT if no reportfile-name is given.
Status and error messages will be written to STDERR if no logfile-name is given.

=head2 Parameters

There are no positional parameters.

The following command-line parameters are supported:

=over 4

=item qza (Q)

Mandatory name of the QZA read-archive to be analyzed.

=item ssu (S)

Mandatory name of the reference-genome "SSU hypervariable-regions" FASTA file,
using genome-IDs as the sequence-IDs,
and with the human-readable genome-name as the "comment" field.

=item reportname (R)

Optinal filename for the output qzaReport-formatted read-counts per genome report;
if no reportname is given, the report will be written to STDOUT.

=item minreads (m)

Minimum required number of reads in a sample worth processing (Default: 100).

=item probe (p) 

Minimum allowed length of 100%-identity in the overlap region between forward (left) and reverse (right) reads
(Default: 25nt)

=item (no)strict

Turn on or off the requirement for "strict" 100%-identity over the full length of the overlap region.
(Default is "nostrict": A 100% match is only required within the "probe" region)

=item logfile

Optional filename to write status and error messages to;
if not given, messages will be written to STDERR.

=item keep (k)

Do not delete the temporary working-directory containing this tool's temporary and intermediate files before exiting.

=item tmpdir

Base directory that the script's temporary working-directory will be created in;
default is the value of the environment-variable $TMPDIR;
if $TMPDIR is undefined, the temporary working-directory will be created in the current working-directory.
The name of the temporary working-directory will be '$TMPDIR/tmp_qza_XXXXXX', where 'XXXXX' is random.
The temporary working-directory will be automatically deleted at the end of the run
unless the '--keep' argument is given.

=item verbose

If specified print detailed progress information to error filehandle.

=item debug

If specified print debugging information to error filehandle.

=item help (h)

Print usage-message and exit.

=back

=cut

use strict;
use warnings;
use Data::Dumper;
use Carp;

use IPC::Run3;
use Time::HiRes qw( time );

use File::Find;
use File::Temp;
use File::Basename qw( fileparse );
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );

use P3Utils;
use SeedUtils;
use gjoseqlib;


my ($rc, @cmd, $cmd, $trouble);

my $opt = P3Utils::script_opts(
    '',
    [ 'qza|Q=s',       'Amplicon QZA archive to be analyzed',            { required => 1  } ],
    [ 'ssu|S=s',       'FASTA of reference SSU hypervariable regions',   { required => 1  } ],
    [],
    [ 'minreads|m=i',  'Minimum required number of reads in a sample (D: 100)', { default => 100 } ],
    [ 'probe|p=i',     'Minimum allowed length of 100% read-overlap region (D: 25nt)',
      { default => 25 } ],
    [ 'strict!',       'Require 100% match over full length of read-overlap (D: nostrict)',
      { default => 0 } ],
    [],
    [ 'reportname|R=s',  'Report-file name (D: STDOUT)',    ],
    [ 'logfile=s',     'Send messages to named logfile (D: STDERR)',           ],
    [],
    [ 'keep|k',        'Keep intermediate and temporary files',               ],
    [ 'tmpdir=s',      'Base-directory that temp directory will be created in (D: $TMPDIR)', { default => $ENV{TMPDIR} } ],
    [ 'verbose|v',     'print progress information',    ],
    [ 'debug',         'print debugging information',   ],
    );
# print($usage->text), exit(0) if $opt->help;

my $qza_file    = $opt->qza();           ### my ($qza_basename) = basename($qza_file, q(.qza));
my $ssu_fasta   = $opt->ssu();
my $reportname  = $opt->reportname();    ### || $qza_basename.q(.qzaReport.tab);

my $minreads    = $opt->minreads();
my $probe       = $opt->probe()  ? q(--probe=).$opt->probe() : q();
my $strict      = $opt->strict() ? q(--strict) : q();

my $keep        = $opt->keep();
my $logfile     = $opt->logfile();
my $tmpdir_base = $opt->tmpdir();
my $verbose     = $opt->verbose() || $ENV{VERBOSE};
my $debug       = $opt->debug()   || $ENV{DEBUG};

my $err_fh = \*STDERR;
if ($logfile) {
    open($err_fh, q(>), $logfile) or die("Could not write-open logfile '$logfile'");
}

if ($debug) {
    $keep = 1;
    $verbose = 1;
 }


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Set up temporary-file directory...
#-----------------------------------------------------------------------
my $tmpdirObj =  File::Temp->newdir( TEMPLATE => q(tmp_qza_XXXXX),
				     DIR => $tmpdir_base,
				     CLEANUP => !$keep
    );
my $tmpdir = $tmpdirObj->dirname();
print STDERR "tmpdir = '$tmpdir'\n\n" if $ENV{DEBUG};

my $tmpdir_tmp = "$tmpdir/tmp";
mkdir($tmpdir_tmp) or die "Could not create temporary working-directory tmpdir_tmp='$tmpdir_tmp'";


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Index and unpack the QZA 'MANIFEST' and Read files...
#-----------------------------------------------------------------------
print $err_fh qq(Indexing and unpacking qza_file='$qza_file'\n) if $verbose;
my $zip = Archive::Zip->new($qza_file);
my @qza_files = $zip->memberNames();;
print $err_fh (qq(INFO:\tQZA archive contents ---\n),
	       Dumper(\@qza_files),
	       qq(\n)
) if $debug;

my ($qza_basedir)   = ($qza_files[0] =~ m{^([^/]+)/});
@qza_files = map { s/$qza_basedir/$tmpdir/; $_ } @qza_files;
print $err_fh (qq(INFO:\tQZA files after mapping to tmpdir ---\n),
	       Dumper(\@qza_files)
) if $debug;

$rc = $zip->extractTree($qza_basedir, $tmpdir);
die("Extraction failed") unless ($rc eq AZ_OK);

my $manifest_file = "$tmpdir/data/MANIFEST";
die "ERROR: Manifest file not found" unless (-s $manifest_file);
my @manifest = &SeedUtils::file_read($manifest_file);
print STDERR (qq(\nMANIFEST:\n), @manifest, qq(\n)) if $verbose;
chomp @manifest;
shift @manifest;   #...drop header-line

my %samples;
my $tmp_data_dir = qq($tmpdir/data);
use constant FORWARD => 0;
use constant REVERSE => 1;
my $line_number = 1;
foreach my $line (@manifest) {
    my ($sample_id, $sample_readfile, $direction) = split /,/, $line;
    if ($direction eq q(forward)) {
	$samples{$sample_id}->[FORWARD] = $sample_readfile;
    }
    elsif ($direction eq q(reverse)) {
	$samples{$sample_id}->[REVERSE] = $sample_readfile;
    }
    else {
	die("ERROR:\t Invalid direction='$direction' at MANIFEST line ", ++$line_number, qq(\n), $line, qq(\n));
    }
}
print $err_fh (q(-) x 72, qq(\n\n)) if $verbose;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Process samples...
#-----------------------------------------------------------------------
my $report_fh = \*STDOUT;
if ($reportname) {
    open($report_fh, q(>), $reportname)
	or die("ERROR: Could not write-open report-file '$reportname'");
}
print $report_fh (join(qq(\t), qw(sample_id  repgen_id  repgen_name  count)), qq(\n));

foreach my $sample_id (sort keys %samples) {
    $trouble = 0;
    print $err_fh "Processing sample '$sample_id'\n" if $verbose;
    
    my $fwd_fastq = qq($tmp_data_dir/).$samples{$sample_id}->[FORWARD];
    my $rev_fastq = qq($tmp_data_dir/).$samples{$sample_id}->[REVERSE];
    unless ((-s $fwd_fastq > 4096) && (-s $rev_fastq > 4096)) {
	print $err_fh ("WARNING:\tFor sample=$sample_id, FWD and/or REV reads missing or too small --- skipping\n",
		       q(-) x 72,
		       qq(\n)
	    );
	next;
    }
    else {
	print $err_fh (qq(Read FASTQ-file sizes:\n),
		       (-s $fwd_fastq).qq(:\t).$fwd_fastq, qq(\n),
		       (-s $rev_fastq).qq(:\t).$rev_fastq, qq(\n\n)
	    ) if $debug;
    }
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Convert forward (left) FASTQ to FASTA...
#-----------------------------------------------------------------------
    print $err_fh "Converting forward reads to FASTA\n" if $verbose;
    my ($fwd_basename)  = fileparse($fwd_fastq, qr/\.fastq(\.gz)?/);
    print $err_fh qq(Could not parse fwd_fastq='$fwd_fastq') and die unless $fwd_basename;
    my $fwd_reads = qq($tmpdir_tmp/$fwd_basename.fna);
    #die Dumper($fwd_fastq, $fwd_basename, $fwd_reads);
    
    print $err_fh qq(Converting '$fwd_fastq' to '$fwd_reads'\n) if $verbose;
    @cmd = (q(seqtk), q(seq), q(-a), $fwd_fastq); # input.fastq > output.fasta
    $cmd = join(q( ), (@cmd, q(>), $fwd_reads));
    $rc  = &run_safe( \@cmd, \undef, $fwd_reads, $err_fh);
    print $err_fh "\n" if $verbose;
    
    if ((-s $fwd_reads) == 0) {
	$trouble = 1;
	print $err_fh "WARNING:\tForward read-FASTA has zero size\n";
    }
    else {
	my $num_reads = map { $_->[0] } &gjoseqlib::read_fasta($fwd_reads);
	if ($num_reads < $minreads) {
	    $trouble = 1;
	    print $err_fh "WARNING:\tForward read-FASTA only contains $num_reads reads\n";
	}
    }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Convert reverse (right) FASTQ to FASTA...
#----------------------------------------------------------------------
    print $err_fh "Converting reverse reads to FASTA\n" if $verbose;
    my ($rev_basename)  = fileparse($rev_fastq, qr/\.fastq(\.gz)?/);
    print $err_fh qq(Could not parse rev_fastq='$rev_fastq') and die unless $rev_basename;
    my $rev_reads = qq($tmpdir_tmp/$rev_basename.fna);
    #die Dumper($rev_fastq, $rev_basename, $rev_reads);
    
    print $err_fh qq(Converting '$rev_fastq' to '$rev_reads'\n) if $verbose;
    @cmd = (q(seqtk), q(seq), q(-a), $rev_fastq); # input.fastq > output.fasta
    $cmd = join(q( ), (@cmd, q(>), $rev_reads));
    $rc  = &run_safe( \@cmd, \undef, $rev_reads, $err_fh);
    print $err_fh "\n" if $verbose;
    
    if ((-s $rev_reads) == 0) {
	$trouble = 1;
	print $err_fh "WARNING:\tReverse read-FASTA has zero size\n";
    }
    else {
	my $num_reads = map { $_->[0] } &gjoseqlib::read_fasta($rev_reads);
	if ($num_reads < $minreads) {
	    $trouble = 1;
	    print $err_fh "WARNING:\tReverse read-FASTA only contains $num_reads reads\n";
	}
    }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Skip sample if trouble    
#-----------------------------------------------------------------------
    if ($trouble) {
	print $err_fh ("INFO:\tSkipping sample_id='$sample_id' because missing or too small\n",
		       q(-) x 72,
		       qq(\n)
	    );
	next;
    }


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Validate required input data...
#-----------------------------------------------------------------------
    my $t0 = time();
    
    -s $fwd_reads
	or print STDERR "Forward reads not found.\n", $opt->help()
	and exit;
    
    -s $rev_reads
	or print STDERR "Reverse reads not found.\n", $opt->help()
	and exit;
    
    my $blast_nsq = $ssu_fasta.q(.nsq);
    -s $blast_nsq
	or print STDERR "SSU rRNA blast database not found.\n", $opt->help()
	and exit;
    my $basename = "$tmpdir_tmp/rRNA_analysis";
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Join overlapping reads...
#-----------------------------------------------------------------------
    print STDERR "\n";
    print STDERR "Running overlap_reads ...\n";
    #system( qq(overlap_reads -s '$fwd_reads' '$rev_reads' > '$basename.joined') );  # -s stringent overlap
    
    @cmd = (q(overlap_reads), ($strict ? $strict : ()), ($probe ? $probe : ()), qq($fwd_reads), qq($rev_reads));
    $cmd = join(q( ), (@cmd, q(>), qq('$basename.joined')));
    print STDERR qq(Running:\t$cmd\n);
    $rc = &run_safe(\@cmd, \undef, qq($basename.joined), $err_fh);
    
    my $t1 = time();
    printf STDERR "... %.3f seconds\n\n", $t1-$t0;
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Make reads nonredundant...
#-----------------------------------------------------------------------
    print STDERR "Running make_nonredundant ...\n";
    #system( qq(make_nonredundant < '$basename.joined' > '$basename.unique') );
    @cmd = ( q(make_nonredundant), q(--input), qq($basename.joined), q(--output), qq($basename.unique));
    
    #$cmd = join(q( ), (@cmd, q(<), qq('$basename.joined'), q(>), qq('$basename.unique')));
    $cmd = join(q( ), @cmd);
    
    print STDERR (qq(RUNNING:\t$cmd\n));
    # run( \@cmd, q(<), qq($basename.joined), q(>), qq($basename.unique) )
    #     or die(qq(ERROR --- Command failed:\t$cmd\n));
    $rc = &run_safe(\@cmd, \undef, \undef, $err_fh);
    
    my $t2 = time();
    printf STDERR "... %.3f seconds\n\n", $t2-$t1;
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...BLASTN of joined reads against database...
#-----------------------------------------------------------------------
    print STDERR "Running blastn ...\n";
    #system( qq(blastn -query '$basename.unique' -db '$ssu_fasta' -evalue 1e-40 -perc_identity 40 -qcov_hsp_perc 80 -penalty -1 -reward 1 -gapopen 2 -gapextend 1 -num_descriptions 5 -num_alignments 5 -num_threads 8 > '$basename.blastn') );
    @cmd = (q(blastn),
	    q(-query),     qq($basename.unique),
	    q(-db),        $ssu_fasta,
	    q(-evalue),         1e-40,
	    q(-perc_identity),     40,
	    q(-qcov_hsp_perc),     80,
	    q(-penalty),           -1,
	    q(-reward),             1,
	    q(-gapopen),            2,
	    q(-gapextend),          1,
	    q(-num_descriptions),   5,
	    q(-num_alignments),     5,
	    q(-num_threads),       32,
	);
    $cmd = join(q( ), (@cmd, q(>), qq('$basename.blastn')));
    print STDERR (qq(RUNNING:\t$cmd\n));
    $rc = &run_safe( \@cmd, \undef, qq($basename.blastn), $err_fh);
    my $t3 = time();
    printf STDERR "... %.3f seconds\n\n", $t3-$t2;
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Compile organism read-counts...
#-----------------------------------------------------------------------
    print STDERR "Running compile_org_counts ...\n";
    #system( qq(compile_org_counts < '$basename.blastn' > '$basename.org_counts') );
    @cmd = qw( compile_org_counts );
    $cmd = join(q( ), (@cmd, q(<), qq('$basename.blastn'), q(>), qq('$basename.org_counts')));
    print STDERR (qq(RUNNING:\t$cmd\n));
    $rc = &run_safe( \@cmd, qq($basename.blastn), qq($basename.org_counts), $err_fh);
    my $t4 = time();
    printf STDERR "... %.3f seconds\n\n", $t4-$t3;
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...Postprocess org-counts into 'qzaReport' format...
#-----------------------------------------------------------------------
    print $report_fh map { chomp;
			   m/^\s+(\S+)\s+(\S+)\s+(.*)$/
			       ? join(qq(\t), ($sample_id, $2, $3, $1)).qq(\n)
			       : ()
    } &SeedUtils::file_read("$basename.org_counts");
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#...End of main loop...
#-----------------------------------------------------------------------
    print $err_fh (q(-) x 72, qq(\n)) if $verbose;
}

print $err_fh "\nTemporary files are in directory '$tmpdir/tmp/'\n" if $keep;
exit(0);


sub run_safe {
#   print STDERR Dumper(\@_);
    my ( $args, $in_fh, $out_fh, $err_fh ) = @_;
#   print $err_fh Dumper($in_fh, $out_fh, $err_fh);
    
    if (my $rc = run3( $args, $in_fh, $out_fh, $err_fh )) {
	return $rc;
    }
    else {
	if ($? == -1) {
	    print $err_fh "ERROR:\tfailed to execute cmd, \$\!=$!\n";
	    confess q(aborted);
	}
	elsif ($? & 127) {
	    print $err_fh ("child died with signal %d, %s coredump\n",
			   (($? & 127),  ($? & 128) ? 'with' : 'without'),
		);
	    confess q(aborted);
	}
    }
}
