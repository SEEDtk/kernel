use strict;
use FIG_Config;
use ScriptUtils;
use Stats;
use GPUtils;
use Math::Round;

$| = 1;
my $stats = Stats->new();
# Start by getting the significant roles.
open(my $ih, "<important_roles.tbl") || die "Could not open important_roles: $!";
my $line = <$ih>;
my %sigRoles = map { ($_ =~ /^(\S+)/) => 1 } <$ih>;
my $totRoles = scalar keys %sigRoles;
print STDERR "$totRoles significant roles found in file.\n";
# Open up GenomePackages and ModPackages.
print STDERR "Loading GenomePackages.\n";
my $gHash = GPUtils::get_all('GenomePackages');
print STDERR "Loading ModPackages.\n";
my $mHash = GPUtils::get_all('ModPackages');
# Write the header.
print join("\t", qw(genome name evalConG evalConF evalGcomplt evalGcontam scikitG scikitF checkMcomplt checkMcontam sigRoles) ) . "\n";
# We read the quality rows for both versions, then propose an alternate score for
# the ModPackages SciKit.
print STDERR "Looping through genomes.\n";
for my $genome (sort keys %$gHash) {
    my $mDir = $mHash->{$genome};
    my $gDir = $gHash->{$genome};
    if (! $mDir) {
        print STDERR "Missing ModPackages directory for $genome.\n";
        $stats->Add(missingGenome => 1);
    } elsif (! -s "$mDir/EvalBySciKit/evaluate.out") {
        $stats->Add(missingEvalCon => 1);
    } elsif (! open(my $qh, "<$mDir/quality.tbl")) {
        print STDERR "Could not open ModPackages quality file for $genome: $!\n";
        $stats->Add(missingMQuality => 1);
    } else {
        my @flds = ScriptUtils::get_line($qh);
        my @cols = ($genome, $flds[2], @flds[11..14]);
        close $qh; undef $qh;
        if (! open($qh, "<$gDir/quality.tbl")) {
            print STDERR "Could not open GenomePackages quality file for $genome: $!\n";
            $stats->Add(missingGQuality => 1);
        } else {
            @flds = ScriptUtils::get_line($qh);
            push @cols, @flds[11..14];
            close $qh; undef $qh;
            if (! open(my $ih, "<$mDir/EvalBySciKit/evaluate.out")) {
                print STDERR "Could not open ModPackages role analysis: $!\n";
                $stats->Add(openFail => 1);
            } else {
                print STDERR "Processing $genome.\n";
                my ($count, $total) = (0, 0);
                while (! eof $ih) {
                    my ($role, $predicted, $actual) = ScriptUtils::get_line($ih);
                    if (! $role) {
                        $stats->Add(badRole => 1);
                    } elsif (defined $predicted && defined $actual && $sigRoles{$role}) {
                        $stats->Add(importantRole => 1);
                        $total++;
                        if ($predicted == $actual) {
                            $count++;
                            $stats->Add(goodRole => 1);
                        }
                    }
                }
                my $result = 0;
                if ($total > 0) {
                    $result = Math::Round::nearest(0.01, $count * 100 / $total);
                }
                push @cols, $result;
                print join("\t", @cols) . "\n";
            }
        }
    }
}
print STDERR "All done.\n" . $stats->Show();