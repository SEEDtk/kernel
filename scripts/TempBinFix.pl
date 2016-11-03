use strict;
use FIG_Config;
use ScriptUtils;
use GenomeTypeObject;

$| = 1;
my $opt = ScriptUtils::Opts('directory');
my ($inDir) = @ARGV;
opendir(my $dh, "$inDir") || die "Could not open input directory $inDir: $!";
my @samples = grep { -s "$inDir/$_/bins.report.txt" } readdir $dh;
for my $sample (@samples) {
    # Read in the whole bin file. Remember the bin headers.
    my $reportFile = "$inDir/$sample/bins.report.txt";
    open(my $ih, "<$reportFile") || die "Could not open report file for $sample: $!";
    my (@reportFile, %headers);
    while (! eof $ih) {
        my $line = <$ih>;
        if ($line =~ /^BIN\s+\d+\s+\(from\s+([^,]+)/) {
            $headers{$1} = scalar(@reportFile);
        }
        push @reportFile, $line;
    }
    close $ih;
    # Loop through the GTO files.
    my $i = 1;
    while (-s "$inDir/$sample/bin$i.gto") {
        my $gto = GenomeTypeObject->create_from_file("$inDir/$sample/bin$i.gto");
        my $contigList = $gto->{contigs};
        # Find a contig that matches a header.
        my $found;
        while (! $found && scalar @$contigList) {
            my $contig = pop @$contigList;
            my $contigID = $contig->{id};
            if (exists $headers{$contigID}) {
                my $idx = $headers{$contigID};
                if ($reportFile[$idx] =~ /^BIN\s+\d+(.+)/) {
                    $reportFile[$idx] = "BIN $i $1\n";
                    $found = 1;
                }
            }
        }
        # Move to the next bin.
        $i++;
    }
    # Write the corrected report file.
    open(my $oh, ">$reportFile") || die "Could not rewrite report file: $!";
    print $oh @reportFile;
    close $oh;
}

