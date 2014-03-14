#! /usr/bin/perl

use strict;

my @sample=();
open (in, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
my $line=<in>;
close in;
my @a=split(/\t/, $line);
for my $i (0..$#a) {
    if ($a[$i] =~ /_class$/) {
	my $name=$a[$i];
	$name =~ s/_class//;
	my $j=$i+1;
	my $k=$i+2;
	my $l=$i+3;
	system("cut -f7,4,6,$j,$k,$l $ARGV[0] > temp");
	open (input, "<temp") or die "Can't open temp since $!\n";
	open (output, ">>$name.insertion.bp.bed") or die "Can't open $name.insertion.bp.bed since $!\n";
	my $header=<input>;
	while (my $line=<input>) {
	    chomp($line);
	    my @b=split(/\t/, $line);
	    if (($b[4] ne "0")||($b[5] ne "0")) {
		my @c=split(/\:/, $b[2]);
		my @d=split(/\./, $c[1]);
		if ($d[0] > $d[1]) {
		    my $temp=$d[0];
		    $d[0]=$d[1];
		    $d[1]=$temp;
		}
		my $lower=$d[0];
		my $upper=$d[1];
		print output "$c[0]\t$lower\t$upper\t$b[0]\t$b[1]\t$b[3]\t$b[4]\t$b[5]\n";
	    }
	}
	close input;
	close output;
	system("rm temp");
    }
}
	
    
