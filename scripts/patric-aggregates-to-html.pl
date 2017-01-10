#!/usr/bin/perl
use strict;
use URI::Escape;
use SeedUtils;
use ScriptUtils;
use Shrub;
use Shrub::Roles;

## Reads standard input.
## Input is from p3-format-results.
my $patric = "https://www.beta.patricbrc.org/";

my $opt = ScriptUtils::Opts('', ScriptUtils::ih_options(), Shrub::script_options());
my $hdg=1;
my $f1;
my $f2;
my $count;
my $genome;
my $html = "";

my $shrub = Shrub->new_for_script($opt);
my $ih = ScriptUtils::IH($opt->input);
while (<$ih>) {
    if ($_ =~ '////') {
        $hdg=1;
    } elsif ($_ =~ '//$') {
         next;
    } elsif ($_ =~ '^###') {
        my ($hash, $genome) = split("\t", $_);
        #print  "<H3>$genome </H3>";
        print  "<table border=\"1\">\n";
        print "<th colspan=3><h3>$genome<h3></th><th>Subsystems</th><th>Reactions</th>\n";
        print $html;
        print "</table><br><br><br>\n\n";
        $html="";
    }else {
        if ($hdg) {
              chomp($_);
              ($f1, $f2, $count) = split("\t",$_);
              print  "<H2> $f1 and $f2 <br>occur together $count times</H2>";
#              print "(<span style=\"color: blue; font-weight: 300;\">fig|nnn.n.peg.n</span> = go to Patric feature page for this peg)";
              print "<br>(<span style=\"color: blue; font-weight: 300;\">&#9400;</span> = go to compare regions for this peg)";
              $hdg=0;
        } else {
            $html .=  "<tr>\n";
            chomp $_;
            my ($id, $fam, $func) = split("\t", $_);
            my ($ss, $reactions) = ('', '');
            if ($func && ! ($func =~ /hypothetical protein/)) {
                $ss = &ss_of_func($func, $shrub);
                $reactions = &reactions_of_func($func, $shrub);
            }

            my $escId = uri_escape($id);
            my $link = $patric."view/Feature/".$escId;
            my $crlink = "http://p3.theseed.org/qa/compare_regions/$escId";
            $html .=  "<td><A HREF=\"".$link."\" target=\_blank >".$id."</A>&nbsp &nbsp";
            $html .= "<A HREF=\"".$crlink."\" target=\_blank style=\"font-size: 100%; font-weight: 300; color: blue;\">&#9400;</A></td>\n";
            my $color = "color:blue";
            if ($fam eq $f1 || $fam eq $f2) {$color="color:red";}
            my $famlink = $patric."view/FeatureList/?eq(plfam_id,$fam)#view_tab=features";
            $html .=  "<td><A HREF=\"".$famlink."\" target=\_blank style=\"$color\">".$fam."</A></td>\n";
            $html .=  "<td>$func</td>\n";

            $html .=  "<td>$ss</td>\n";
            $html .=  "<td>$reactions</td>\n";

            $html .=  "</tr>\n";
        }
    }
}

sub ss_of_func {
    my($func, $shrub) = @_;

    my @roles = &SeedUtils::roles_of_function($func);
    my %in_ss;
    for my $role (@roles) {
        my $checksum = Shrub::Roles::Checksum($role);
        my @ss = $shrub->GetFlat('Role Subsystem', 'Role(checksum) = ?', [$checksum], 'Subsystem(name)');
        for my $ss (@ss) {
            $in_ss{$ss} = 1;
        }
    }
    return join("\n",map { "<p>$_</p>" }sort keys %in_ss);
}

sub reactions_of_func {
    my($func, $shrub) = @_;

    my @roles = &SeedUtils::roles_of_function($func);
    my %has_react;
    for my $role (@roles) {
        my $checksum = Shrub::Roles::Checksum($role);
        my @react = $shrub->GetFlat('Role Complex2Reaction', 'Role(checksum) = ?', [$checksum], 'Complex2Reaction(to-link)');
        for my $react (@react) {
            if (! $has_react{$react}) {
                $has_react{$react} = $shrub->reaction_formula($react);
            }
        }
    }
    return join("\n",map { "<P>$_</P>" } sort values %has_react);
}
