use strict;
use FIG_Config;
use ScriptUtils;
use Stats;

open(my $ih, "<Cancer/dose_response/GDSC_IC50.tbl") || die "Could not open GDSC input: $!";
open(my $oh, ">Cancer/dose_response/GDSC_IC50") || die "Could not open GDSC output: $!";
my $line = <$ih>;
print $oh join("\t", qw(drug cell_line ic50)) . "\n";
while (! eof $ih) {
    my @cols = ScriptUtils::get_line($ih);
    my $drug = "GDSC.$cols[4]";
    my $cl = "GDSC.$cols[3]";
    my $ic50 = - $cols[9] / log(10);
    print $oh join("\t", $drug, $cl, $ic50) . "\n";
}
close $ih; undef $ih;
close $oh; undef $oh;
open($ih, "<Cancer/dose_response/NCI60_IC50.tbl") || die "Could not open NCI60 input: $!";
open($oh, ">Cancer/dose_response/NCI60_IC50") || die "Could not open NCI60 output: $!";
print $oh join("\t", qw(drug cell_line ic50)) . "\n";
$line = <$ih>;
while (! eof $ih) {
    my @cols = ScriptUtils::get_line($ih);
    my $drug = "NSC.$cols[0]";
    my $cl = "NCI60.$cols[4]";
    my $ic50 = - $cols[7];
    print $oh join("\t", $drug, $cl, $ic50) . "\n";
}
