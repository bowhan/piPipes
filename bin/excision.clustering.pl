#! /usr/bin/perl

use strict;

my %position=();
my %names=();
open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    my @b=split(/\#/, $a[8]);

    if (defined $position{$b[0]}) {
	my @c=split(/\:/, $position{$b[0]});
	if (($c[0] eq $a[0])&&($a[1] < $c[1])) {
	    $position{$b[0]} =~ s/$c[1]/$a[1]/;
	}
	if (($c[0] eq $a[0])&&($a[2] > $c[2])) {
	    $position{$b[0]} =~ s/$c[2]/$a[2]/;
	}
	my $transposon=$a[3];
	if ($names{$b[0]} !~ /$transposon/) {$names{$b[0]}=$names{$b[0]}.",$transposon";}
    }
    else {
	$position{$b[0]}="$a[0]\:$a[1]\:$a[2]";
	$names{$b[0]}=$a[3];
    }
}
close input;

open (output, ">>temp_for_sort") or die "Can't open temp_for_sort since $!\n";
while ((my $key, my $value) = each (%position)) {
    my @z=split(/\:/, $value);
    print output "$z[0]\t$z[1]\t$z[2]\t$names{$key}\n";
}
close output;

system("sort +0 -1 +1n -2 +2n -3 temp_for_sort > sorted");
system("uniq -c sorted > $ARGV[1]");
system("rm sorted temp_for_sort");
