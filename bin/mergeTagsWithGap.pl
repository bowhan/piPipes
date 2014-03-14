#!/share/bin/perl
#chr2L   114333  114409  FBgn0003055_P-element   harwich,1;      +       antisense
#chr2L   114443  114567  FBgn0003055_P-element   harwich,3;      +       antisense
#chr2L   114636  114712  FBgn0003055_P-element   harwich,1;      -       antisense
#chr2L   131640  131929  FBgn0003055_P-element   harwich,42;     +       sense
#chr2L   131948  132274  FBgn0003055_P-element   harwich,18;     -       sense
#chr2L   132027  132103  FBgn0003055_P-element   harwich,1;      -       antisense

use warnings;
use strict;
use List::Util qw(max min);

if(scalar(@ARGV)<2 || grep {/^-h/} @ARGV)
{
	die "
usage: mergeOverlapBed4.pl inputFile
Expects BED input with at least 4 fields.  For each {chr,name} pair,
merges overlapping ranges and prints out sorted BED4 to stdout.
inputFile can be - or stdin to read from stdin.
";
}

my $input=shift @ARGV;
my $maxgap=shift @ARGV;
grep {s/^stdin$/-/i} $input;

my %item2coords;
open IN,$input;
while (<IN>)
{
	chomp;
	my ($chrom,$start,$end,$transposonName,$count,$strand,$transposonStrand)=split/\t/;
	push @{$item2coords{"$chrom;$transposonName;$transposonStrand"}},[$start,$end,$count,$strand]; 
}
close IN;

my @results;
foreach my $item (keys %item2coords)
{
	my @sortedCoords=sort{ $a->[0]<=>$b->[0] } @{$item2coords{$item}};
	my ($chrom,$tName,$tStrand)=split(/;/,$item);
	my ($mergeStart,$mergeEnd,$mergeCounts,$mergeStrand)=@{shift @sortedCoords};
	my %sampleCounts=();
	my ($breakStart,$breakEnd)=0;
	foreach my $sa (split/;/,$mergeCounts)
	{
		my ($s,$c)=split/,/,$sa;
		$sampleCounts{$s}{$mergeStrand}=$c;
	}
	foreach my $rangeRef (@sortedCoords) 
	{
    		my ($rangeStart,$rangeEnd,$rangeCounts,$rangeStrand)=@{$rangeRef};
		if($mergeStrand=~/\Q$rangeStrand\E$/)
		{
			if($rangeStart>=$mergeEnd+$maxgap)
			{
				$mergeCounts="";
				foreach my $s (keys %sampleCounts)
				{
					$mergeCounts.=$s;
					$mergeCounts.=",".$_.$sampleCounts{$s}{$_} foreach keys %{$sampleCounts{$s}};
					$mergeCounts.=";";
				}
				if($mergeStrand eq "+")
				{
					$breakStart=$mergeEnd;
					$breakEnd=$mergeEnd+$maxgap;
				}
				if($mergeStrand eq "-")
				{
					$breakStart=$mergeStart-$maxgap;
					$breakEnd=$mergeStart;
				}
				push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$mergeCounts,$mergeStrand,$tStrand,"$chrom:$breakStart.$breakEnd"];
				($mergeStart,$mergeEnd,$mergeStrand)=($rangeStart,$rangeEnd,$rangeStrand);
				%sampleCounts=();
				foreach my $sa (split/;/,$rangeCounts)
				{
					my ($s,$c)=split/,/,$sa;
					$sampleCounts{$s}{$rangeStrand}=$c;
				}
			}
			else
			{
				$mergeEnd=max($rangeEnd,$mergeEnd);
				foreach my $sa (split/;/,$rangeCounts)
				{
					my ($s,$c)=split/,/,$sa;
					$sampleCounts{$s}{$rangeStrand}+=$c;
				}
			}
		}
		elsif($rangeStrand eq "+")
		{
			$mergeCounts="";
			foreach my $s (keys %sampleCounts)
			{
				$mergeCounts.=$s;
				$mergeCounts.=",".$_.$sampleCounts{$s}{$_} foreach keys %{$sampleCounts{$s}};
				$mergeCounts.=";";
			}
			if($mergeStrand eq "+")
			{
				$breakStart=$mergeEnd;
				$breakEnd=$mergeEnd+$maxgap;
			}
			if($mergeStrand eq "-")
			{
				$breakStart=$mergeStart-$maxgap;
				$breakEnd=$mergeStart;
			}
			push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$mergeCounts,$mergeStrand,$tStrand,"$chrom:$breakStart.$breakEnd"];
			($mergeStart,$mergeEnd,$mergeStrand)=($rangeStart,$rangeEnd,$rangeStrand);
			%sampleCounts=();
			foreach my $sa (split/;/,$rangeCounts)
			{
				my ($s,$c)=split/,/,$sa;
				$sampleCounts{$s}{$rangeStrand}=$c;
			}
		}
		else
		{
			if($rangeStart>=$mergeEnd+$maxgap*2)
			{
				$mergeCounts="";
				foreach my $s (keys %sampleCounts)
				{
					$mergeCounts.=$s;
					$mergeCounts.=",".$_.$sampleCounts{$s}{$_} foreach keys %{$sampleCounts{$s}};
					$mergeCounts.=";";
				}
				if($mergeStrand eq "+")
				{
					$breakStart=$mergeEnd;
					$breakEnd=$mergeEnd+$maxgap;
				}
				if($mergeStrand eq "-")
				{
					$breakStart=$mergeStart-$maxgap;
					$breakEnd=$mergeStart;
				}
				push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$mergeCounts,$mergeStrand,$tStrand,"$chrom:$breakStart.$breakEnd"];
				($mergeStart,$mergeEnd,$mergeStrand)=($rangeStart,$rangeEnd,$rangeStrand);
				%sampleCounts=();
				foreach my $sa (split/;/,$rangeCounts)
				{
					my ($s,$c)=split/,/,$sa;
					$sampleCounts{$s}{$rangeStrand}=$c;
				}
			}
			else
			{
				$breakStart=$mergeEnd;
				$mergeEnd=max($rangeEnd,$mergeEnd);
				$breakEnd=$rangeStart;
				foreach my $sa (split/;/,$rangeCounts)
				{
					my ($s,$c)=split/,/,$sa;
					$sampleCounts{$s}{$rangeStrand}+=$c;
				}
				$mergeStrand.=$rangeStrand;
			}			
		}
	}
	$mergeCounts="";
	foreach my $s (keys %sampleCounts)
	{
		$mergeCounts.=$s;
		$mergeCounts.=",".$_.$sampleCounts{$s}{$_} foreach keys %{$sampleCounts{$s}};
		$mergeCounts.=";";
	}
	if($mergeStrand eq "+")
	{
		$breakStart=$mergeEnd;
		$breakEnd=$mergeEnd+$maxgap;
	}
	if($mergeStrand eq "-")
	{
		$breakStart=$mergeStart-$maxgap;
		$breakEnd=$mergeStart;
	}
	push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$mergeCounts,$mergeStrand,$tStrand,"$chrom:$breakStart.$breakEnd"] if $mergeEnd;
}

sub bed4Cmp
{
  # For sorting by chrom, chromStart, and names -- reverse order for names
	return $a->[0] cmp $b->[0] ||
	$a->[1] <=> $b->[1] ||
	$b->[3] cmp $a->[3];
}

foreach my $r (sort bed4Cmp @results)
{
	print join("\t",@{$r}),"\n";
}
