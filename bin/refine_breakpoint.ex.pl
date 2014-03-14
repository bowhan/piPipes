#! /usr/bin/perl

use strict;

my @files=<*.excision.cluster.*>;
foreach my $file (@files) {
    if (($file !~ /sfcp/)&&($file !~ /refsup/)) {
	my $sfcp=$file.".sfcp";
	my $title=$file.".refined.bp";
	
	open (input, "<$file") or die "Can't open $file since $!\n";
	open (input1, "<$sfcp") or die "Can't open $sfcp since $!\n";
	open (output, ">>$title") or die "Can't open $title since $!\n";
	print output "Chr\tStart\tEnd\tTransposonName\t5\'_Junction\t3\'_Junction\n";
	while (my $line=<input>) {
	    chomp($line);
	    my @a=split(/\s+/, $line);
	    my $line1=<input1>;
	    chomp($line1);
	    my @b=split(/\t/, $line1);
	    my @pos=split(/\;/, $b[4]);
	    my $plusnext=""; my $minusnext="";
	    my $plusover=0; my $minusover=0;
	    my $lpcoor=""; my $lmcoor=""; my $rpcoor=""; my $rmcoor="";
	    my $lp=0; my $lm=0; my $rp=0; my $rm=0;
	    my %plus=(); my %minus=();
	    foreach my $site (@pos) {
		my @x=split(/\:/, $site);
		my @y=split(/\(/, $x[0]);
		chop($y[1]);
		if (($y[0] =~ /\-/)&&($y[1] eq "+")&&($x[1] >= $plusover)) {
		    if ($plusover >= 2) {$plusnext="$lpcoor\-$rpcoor\:$plusover";}
		    $plusover=$x[1]; 
		    my @z=split(/\-/, $y[0]);
		    $lpcoor=$z[0]; $lp=$x[1];
		    $rpcoor=$z[1]; $rp=$x[1];
		}
		elsif (($y[0] =~ /\-/)&&($y[1] eq "-")&&($x[1] >= $minusover)) {
		    if ($minusover >= 2) {$minusnext="$lmcoor\-$rmcoor\:$minusover";}
		    $minusover=$x[1];
		    my @z=split(/\-/, $y[0]);
		    $lmcoor=$z[0]; $lm=$x[1];
		    $rmcoor=$z[1]; $rm=$x[1];
		}
		elsif (($y[0] !~ /\-/)&&($y[1] eq "+")) {
		    $plus{$y[0]}=$x[1];
		}
		elsif (($y[0] !~ /\-/)&&($y[1] eq "-")) {
		    $minus{$y[0]}=$x[1];
		}
	    }
	    
	    if (($plusnext ne "")&&($minusover == 0)) {
		my @m=split(/\:/, $plusnext);
		if (($m[1] >= 2)&&($m[1] == $plusover)) {
		    my $count1=$m[1]; my $count2=$plusover;
		    foreach my $id (keys %plus) {
			if ($m[0] =~ /$id/) {$count1 += $plus{$id};}
			elsif ($rpcoor == $id) {$count2 += $plus{$id};}
		    }
		    if ($count1 > $count2) {
			my @n=split(/\-/, $m[0]);
			$lpcoor=$n[0]; $lp=$m[1];
			$rpcoor=$n[1]; $rp=$m[1];
		    }
		}
	    }

	    if (($minusnext ne "")&&($plusover == 0)) {
		my @m=split(/\:/, $minusnext);
		if (($m[1] >= 2)&&($m[1] == $minusover)) {
		    my $count1=$m[1]; my $count2=$minusover;
		    foreach my $id (keys %minus) {
			if ($m[0] =~ /$id/) {$count1 += $minus{$id};}
			elsif ($lmcoor == $id) {$count2 += $minus{$id};}
		    }
		    if ($count1 > $count2) {
			my @n=split(/\-/, $m[0]);
			$lmcoor=$n[0]; $lm=$m[1];
			$rmcoor=$n[1]; $rm=$m[1];
		    }
		}
	    }			    

	    if (($plusover >= 2)&&($minusover >= 2)&&(($lpcoor-$rpcoor) != ($lmcoor-$rmcoor))) {
		if ($plusnext ne "") {
		    my @m=split(/\:/, $plusnext);
		    my @n=split(/\-/, $m[0]);
		    if ((($n[1]-$n[0]) == ($rmcoor-$lmcoor))&&($m[1] >= 2)) {
			$rpcoor=$n[1];
			$lpcoor=$n[0];
			$plusover=$m[1];
			$lp=$m[1];
			$rp=$m[1];
		    }
		}
		if ($minusnext ne "") {
                    my @m=split(/\:/, $minusnext);
                    my @n=split(/\-/, $m[0]);
                    if ((($n[1]-$n[0]) == ($rpcoor-$lpcoor))&&($m[1] >= 2)) {
			$rmcoor=$n[1];
			$lmcoor=$n[0];
			$minusover=$m[1];
			$lm=$m[1];
			$rm=$m[1];
                    }
		}
	    }

	    my $plusc=0; my $pluscoor="";
	    my $minusc=0; my $minuscoor="";
	    foreach my $id (keys %plus) {
		if ($id eq $rpcoor) {
		    $rp=$plusover+$plus{$id};
		}
		if ($plus{$id} > $plusc) {
		    $plusc=$plus{$id};
		    $pluscoor=$id;
		}
		elsif (($plus{$id} == $plusc)&&(abs($id-$b[2]) < abs($pluscoor-$b[2]))) {
                    $plusc=$plus{$id};
                    $pluscoor=$id;
		}
	    }
	    foreach my $id (keys %minus) {
		if ($id eq $lmcoor) {
		    $lm=$minusover+$minus{$id};
		}
		if ($minus{$id} > $minusc) {
		    $minusc=$minus{$id};
		    $minuscoor=$id;
		}
		elsif (($minus{$id} == $minusc)&&(abs($id-$b[1]) < abs($minuscoor-$b[1]))) {
                    $minusc=$minus{$id};
                    $minuscoor=$id;
		}
	    }
	    if ($plusover < 2) {
		$lpcoor="";
		if ($plusc >= 3) {$rpcoor=$pluscoor; $rp=$plusc;}
		else {$rpcoor="";}
	    }
	    if ($minusover < 2) {
		$rmcoor="";
		if ($minusc >= 3) {$lmcoor=$minuscoor; $lm=$minusc;}
		else {$lmcoor="";}
	    }	    
		
	    my $bp1=""; my $bp2="";
	    if (($lpcoor ne "")&&($lmcoor ne "")) {
		$bp1="$lpcoor\(\+\)\:$lp,$lmcoor\(\-\)\:$lm";
	    }
	    elsif ($lpcoor ne "") {$bp1="$lpcoor\(\+\)\:$lp";}
	    elsif ($lmcoor ne "") {$bp1="$lmcoor\(\-\)\:$lm";}
	    if (($rpcoor ne "")&&($rmcoor ne "")) {
		$bp2="$rpcoor\(\+\)\:$rp,$rmcoor\(\-\)\:$rm";
	    }
	    elsif ($rpcoor ne "") {$bp2="$rpcoor\(\+\)\:$rp";}
	    elsif ($rmcoor ne "") {$bp2="$rmcoor\(\-\)\:$rm";}

	    print output "$a[2]\t$a[3]\t$a[4]\t$a[5]\t$bp1\t$bp2\n";
	}
	
	close input;
	close input1;
	close output;
    }
}
