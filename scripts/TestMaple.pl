use strict;
use P3DataAPI;
use P3Utils;
use GPUtils;
use Stats;
use GenomeTypeObject;

$| = 1;
# Get the command-line options.
my $opt = P3Utils::script_opts('gtoDir', P3Utils::ih_options(),
        ['qual', 'compute standard quality measure'],
        ['recursive', 'search subdirectories for GTOs'],
        );
my $stats = Stats->new();
# Get the working directory.
my ($gtoDir) = @ARGV;
if (! $gtoDir) {
    die "No input directory specified.";
} elsif (! -d $gtoDir) {
    die "Input directory $gtoDir not found or invalid.";
}
# Get the options.
my $qual = $opt->qual;
my $recursive = $opt->recursive;
# Locate the gtos.
print "Finding GTOs.\n";
my @subDirs;
if ($recursive) {
    opendir(my $dh, $gtoDir) || die "Could not open $gtoDir: $!";
    @subDirs = map { "$gtoDir/$_" } grep { substr($_,0,1) ne '.' && -d "$gtoDir/$_" } readdir $dh;
} else {
    push @subDirs, $gtoDir;
}
my %gtos;
for my $subDir (@subDirs) {
    print "Searching $subDir.\n";
    opendir(my $dh, $subDir) || die "Could not open $subDir: $!";
    while (my $file = readdir $dh) {
        if (-s "$subDir/$file" && $file =~ /^(\d+\.\d+)\.gto$/) {
            $gtos{$1} = "$subDir/$file";
        }
    }
}
print scalar(keys %gtos) . " GTOs found in $gtoDir.\n";
# Open the input file.
print "Connecting to input file.\n";
my $ih = P3Utils::ih($opt);
# Read in the header line and find the critical columns.
my @qual;
if ($qual) {
    @qual = ('completeness', 'contamination');
}
my ($headers, $cols) = P3Utils::find_headers($ih, input => 'PATRIC ID', 'Good?', @qual);
my @outHeaders = @$headers;
push @outHeaders, qw(Fids EvalG_group GTO GTO_species);
if ($qual) {
    push @outHeaders, 'Quality';
}
# Open the main output file.
open(my $oh, ">$gtoDir/genomes.tbl") || die "Could not open main output: $!";
# Write the new headers.
P3Utils::print_cols(\@outHeaders, oh => $oh);
# Loop through the input.
my $count = 0;
while (! eof $ih) {
    # Get the next line.
    my $line = <$ih>;
    my @fields = P3Utils::get_fields($line);
    $stats->Add(lineIn => 1);
    $count++;
    # Get the genome ID and goodness flag.
    my ($genomeID, $goodFlag, $complete, $contam) = P3Utils::get_cols(\@fields, $cols);
    # We need to process the seed protein.
    my $evalType = ($goodFlag ? 'good' : 'bad');
    # These will be our extra fields.
    my ($pegs, $qGroup, $species) = (0, '', '');
    # Read the GTO and get its seed protein.
    my $gtoFile = $gtos{$genomeID};
    if (! $gtoFile) {
        print "GTO not found for $genomeID.\n";
        $stats->Add(missingGto => 1);
        $gtoFile = '';
    } else {
        print "Reading GTO for $genomeID ($count).\n";
        my $gto = GenomeTypeObject->create_from_file($gtoFile);
        my $seedProt = GPUtils::get_seed($gto);
        # Verify that the genome is genuinely good.
        if (! $seedProt) {
            $stats->Add(badSeed => 1);
            $fields[$cols->[1]] = '';
            $evalType = 'bad';
        } else {
            $stats->Add(goodSeed => 1);
        }
        my $flist = $gto->{features};
        $pegs = scalar @$flist;
        $qGroup = $gto->{quality}{completeness_group};
        if (! defined $qGroup) {
            $stats->Add(groupMissing => 1);
            $qGroup = '';
        }
        $species = $gto->{ncbi_species};
        if (! defined $species) {
            $stats->Add(speciesMissing => 1);
            $species = '';
        }
    }
    push @fields, $pegs, $qGroup, $gtoFile, $species;
    # Compute the alternate quality measure.
    if ($qual) {
        my $qualType = 'L';
        if ($contam <= 5) {
            if ($complete >= 90) {
                $qualType = 'H';
            } elsif ($complete >= 50) {
                $qualType = 'M';
            }
        }
        $stats->Add("$evalType-$qualType" => 1);
        $stats->Add("genome-$qualType" => 1);
        push @fields, $qualType;
    }
    # Write the output record.
    P3Utils::print_cols(\@fields, oh => $oh);
    $stats->Add(lineOut => 1);
    print "$count genomes processed.\n" if $count % 500 == 0;
}
close $oh;
print "All done.\n" . $stats->Show();