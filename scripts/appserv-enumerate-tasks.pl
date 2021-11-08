use Bio::KBase::AppService::Client;
use Getopt::Long::Descriptive;
use strict;
use Data::Dumper;

my($opt, $usage) = describe_options("%c %o",
                    ["url|u=s", "Service URL"],
                    ["limit=i", "maximum number of tasks to return (approximate)", { default => 1000 }],
                    ["verbose|v", "Show verbose output from service"],
                    ["help|h", "Show this help message"]);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 0);

my $client = Bio::KBase::AppService::Client->new($opt->url);

my $limit = 25;
my $offset = 0;
my $maximum = $opt->limit;
my $done = 0;

my $count = 0;

while (! $done) {
    my $tasks = $client->enumerate_tasks($offset, $limit);
    if (! scalar @$tasks) {
        $done = 1;
    } else {
        $offset += $limit;
        for my $task (@$tasks) {
            print join("\t", $task->{id}, $task->{app}, $task->{workspace}, $task->{status}, $task->{submit_time}, $task->{completed_time}), "\n";
        }
        $count += scalar @$tasks;
        if ($count >= $maximum) {
            $done = 1;
        }
    }
}

