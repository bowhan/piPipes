#! /usr/bin/perl

use strict;

open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
my %leng=();
my %trans=();
my %coordinate=();
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);
    if (defined $leng{$a[9]}) {
	$trans{$a[9]} += $a[2]-$a[1];
    }
    else {
	$leng{$a[9]}=$a[8]-$a[7]-10;
	$trans{$a[9]}=$a[2]-$a[1];
	$coordinate{$a[9]}="$a[6]\:$a[7]\:$a[8]";
    }
}
close input;

open (output, ">>$ARGV[2]") or die "Can't open $ARGV[2] since $!\n";
while ((my $key, my $value) = each (%coordinate)) {
    if ((($leng{$key}-$trans{$key}) <= $ARGV[1])&&(($leng{$key}-$trans{$key}) >= 0)) {
#    if (($leng{$key}-$trans{$key}) <= 500) {
	my @b=split(/\:/, $value);
	print output "$b[0]\t$b[1]\t$b[2]\t$key\n";
    }
}
close output;
