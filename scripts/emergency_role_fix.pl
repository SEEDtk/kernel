use strict;
use FIG_Config;
use ScriptUtils;
use ProtFamRepo;
use Stats;
use GPUtils;
use P3DataAPI;
use Digest::MD5;
use Encode qw(encode_utf8);
use Shrub;

open(my $ih, "<$FIG_Config::global/roles.in.subsystems") || die "Could not open input roles.in.subsystems: $!";
open(my $oh, "<$FIG_Config::data/roles.in.subsystems") || die "Could not open output roles.in.subsystems: $!";
while (! eof $ih) {
    my $line = <$ih>;
    chomp $line;
    my ($id, $checksum, $role) = split /\t/, $line;
    my $newCheck = RoleParse::Checksum($role);
    if ($newCheck ne $checksum) {
        print "Correcting $id: $role.\n";
    }
    print $oh join("\t", $id, $checksum, $role) . "\n";
}
