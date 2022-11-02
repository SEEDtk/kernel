use strict;
use FIG_Config;
use File::Copy::Recursive;
use Getopt::Long::Descriptive;
use P3Utils;

my $opt = P3Utils::script_opts('a b c');
print join("\n", "", @INC, "");

