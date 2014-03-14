#! /usr/bin/perl

use strict;

#my @files=<*.bp.bed.sfcp>;
my @files=<*.clipped.reads.aln>;
for my $file (@files) {

    my $title=$file;
    $title =~ s/clipped.reads.aln/insertion.refined.bp/;
    my $title2=$file;
    $title2 =~ s/clipped.reads.aln/insertion.bp.bed/;

    open (input, "<$file") or die "Can't open $file since $!\n";
    open (input2, "<$title2") or die "Can't open $title2 since $!\n";
    open (output, ">>$title") or die "Can't open $title since $!\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	my @b=split(/\;/, $a[4]);
	my $plusmax="";
	my $minusmax="";
	my $plus=0;
	my $minus=0;
	my $bp="";

	my $line2=<input2>;
	my @z=split(/\t/, $line2);
	my $psup=abs($z[6]);
	my $msup=abs($z[7]);
	my $strand=$z[4];
	my $class=$z[5];

	for my $element (@b) {
	    my @c=split(/\:/, $element);
	    chop($c[0]);
	    my @d=split(/\(/, $c[0]);
	    if (($d[1] eq "+") && ($c[1] > $plus)) {
		$plusmax=$d[0];
		$plus=$c[1];
	    }
	    elsif (($d[1] eq "-") && ($c[1] > $minus)) {
		$minusmax=$d[0];
		$minus=$c[1];
	    }
	}

	if ((($a[0] =~ /^chr/) || ($a[0] =~ /^\d{1,2}$/) || ($a[0] eq "X") || ($a[0] eq "Y")) && ($a[1] > 0)) {
	    $a[1] += 15;
	    $a[2] -= 15;
	    print output "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$strand\t$class\t";
	    if (($minus >= 1)&&($plus >= 1)&&(abs($plusmax-$minusmax) <= 25)) {
#		$psup += $minus;
#		$msup += $plus;
		print output "$plusmax\(\+\)\t$minusmax\(\-\)\t$plus\t$minus\t";
	    }
	    elsif (($plus >= $minus)&&($plus >= 2)&&($plusmax >= $a[1])&&($plusmax <= $a[2])) {
#		$msup += $plus;
		print output "$plusmax\(\+\)\t$plusmax\(\-\)\t$plus\t0\t";
	    }
	    elsif (($minus >= 2)&&($minusmax >= $a[1])&&($minusmax <= $a[2])) {
#		$psup += $minus;
		print output "$minusmax\(\+\)\t$minusmax\(\-\)\t0\t$minus\t";
	    }
	    else {
		my $mid=int(($a[1] + $a[2])/2);
		print output "$mid\(\+\)\t$mid\(\-\)\t0\t0\t";
	    }
	    print output "$psup\t$msup\n";
	}
    }
    close input;
    close input2;
    close output;
    system("uniq $title > temp");
    system("mv temp $title");
}
