use strict;
use FIG_Config;
use GenomeTypeObject;

opendir(my $dh, "Bins_HMP") || die "Could not open Bins_HMP: $!";
my @samples = grep { -s "Bins_HMP/$_/bin1.gto" } readdir $dh;
print STDERR scalar(@samples) . " samples found.\n";
closedir $dh;
print "sample\tsource\tgenome_id\tgenus\tinconsistent\tincomplete\tcontaminated\tother\tgood\timproved\n";
for my $sample (sort @samples) {
    my $source = getSource($sample);
    my $qual = getQual($sample);
    print STDERR "Processing $source sample $sample.\n";
    opendir(my $sh, "Bins_HMP/$sample") || die "Could not open $sample directory: $!";
    my @bins = sort map { "Bins_HMP/$sample/$_" } grep { $_ =~ /^bin\d+\.gto$/ } readdir $sh;
    closedir $sh;
    print STDERR scalar(@bins) . " bins found.\n";
    for my $bin (@bins) {
        print STDERR "Analyzing $bin.\n";
        printBin($sample, $source, $qual, $bin);
    }
}

sub getQual {
    my ($sample) = @_;
    my %retVal;
    open(my $ih, '<', "Bins_HMP/$sample/Eval/index.tbl") || die "Could not open quality file for $sample: $!";
    my $line = <$ih>;
    while (! eof $ih) {
        $line = <$ih>;
        chomp $line;
        my @fields = split /\t/, $line;
        my $inconsistent = ($fields[9] < 87 ? "1" : "");
        my $incomplete = ($fields[10] < 80 ? "1" : "");
        my $contaminated = ($fields[11] > 10 ? "1" : "");
        my $good = ($fields[14] ? "1" : "");
        my $other = ((! $good && ! ($inconsistent || $incomplete || $contaminated)) ? "1" : "");
        $retVal{$fields[1]} = [$inconsistent, $incomplete, $contaminated, $other, $good];
    }
    return \%retVal;
}

sub getSource {
    my ($sample) = @_;
    my $retVal = "unknown";
    if (open(my $ih, '<', "Bins_HMP/$sample/site.tbl")) {
        my $line = <$ih>;
        ($retVal) = ($line =~ /\t(\S+)/);
    }
    return $retVal;
}

sub printBin {
    my ($sample, $source, $qual, $bin) = @_;
    my $gto = GenomeTypeObject->create_from_file($bin);
    my $id = $gto->{id};
    my $q = $qual->{$id};
    if (! $q) {
        print STDERR "Bin $id is invalid-- skipped.\n";
    } else {
        my $genus = $gto->{ncbi_genus};
        my $a = $gto->{analysis_events};
        my @tools = grep { $_->{tool_name} eq "p3x-improve-gto" } @$a;
        my $improved = (scalar(@tools) ? "1" : "");
        print join("\t", $sample, $source, $id, $genus, @$q, $improved) . "\n";
    }
}
