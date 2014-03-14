#! /usr/bin/perl

use strict;

my @files=<*.excision.cluster.*.refined.bp>;
foreach my $file (@files) {
    my $rfsp=$file.".refsup";
    my $count=$file;
    my $title=$file;
    $count =~ s/.refined.bp//;
    $title =~ s/excision/absence/;
    $title =~ s/.cluster.rpmk//;
    $title .= ".summary";
	
    open (input, "<$file") or die "Can't open $file since $!\n";
    open (input1, "<$count") or die "Can't open $count since $!\n";
    open (input2, "<$rfsp") or die "Can't open $rfsp since $!\n";
    open (output, ">>$title") or die "Can't open $title since $!\n";
    my $header=<input>; 
    chomp($header);
    print output "$header\tVariant\tReference\tFrequency\n";
    while (my $line=<input>) {
	chomp($line);
	my @a=split(/\t/, $line);
	my $line1=<input1>;
	my $line2=<input2>;
	chomp($line1);
	chomp($line2);
	my @b=split(/\s+/, $line1);
	my @c=split(/\t/, $line2);

	my $variant=$b[1];
	my @x=split(/\:/, $a[4]);
	my @y=split(/\:/, $a[5]);
	if ($a[4] =~ /\,/) {
	    my @m=split(/\,/, $x[1]);
	    $variant += $m[0]+$x[2];
	}
	else {$variant += $x[1];}
	if ($a[5] =~ /\,/) {
	    my @n=split(/\,/, $y[1]);
	    $variant += $n[0]+$y[2];
	}
	else {$variant += $y[1];}
	my $ratio=sprintf("%.4f", ($variant*2)/($variant*2+$c[6]));
	
	$line =~ s/:\d+//g;
	print output "$line\t$variant\t$c[6]\t$ratio\n";
    }
    
    close input;
    close input1;
    close input2;
    close output;
}
