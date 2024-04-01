#!/usr/bin/perl

open(my $ih, '<', "ss_good.txt") || die "Could not open input file: $!";
my %names;
my $line = <$ih>;
while (! eof $ih) {
	$line = <$ih>;
	chomp $line;
	if ($names{$line}) {
		print "Duplicate: $line\n";
	} else {
		$names{$line} = 1;
	}
}