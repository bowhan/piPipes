#!/share/bin/perl
#chr2L   735929  736005  HWUSI-EAS1533_0002:1:73:4665:12371#0/2  FBgn0000155_roo,-58;FBgn0000155_roo,-8722;      -
use warnings;
use strict;

if(scalar(@ARGV)<1 || grep {/^-h/} @ARGV)
{
	die "
usage: mergeOverlapBed4.pl inputFile
Expects BED input with at least 4 fields.  For each {chr,name} pair,
merges overlapping ranges and prints out sorted BED4 to stdout.
inputFile can be - or stdin to read from stdin.
";
}

my $input=shift @ARGV;
grep {s/^stdin$/-/i} $input;

my %item2coords;
open IN,$input;
while (<IN>)
{
	chomp;
	my ($chrom,$start,$end,$sample,$class,$strand)=split/\t/;
  	die "Sorry, input must have at least 4 fields of BED.\n" if ! $class;
	# random choose one
#	my @loc=$class=~/(.*?),(\+|-)(.*)/;
#	my $transposonStrand=($strand eq $loc[1])?"antisense":"sense";
#	push @{$item2coords{"$chrom;$strand;$loc[0];$transposonStrand"}},[$start,$end,$sample] 

	# norm by class
	my @loc=map { [/(.*?),(\+|-)(.*)/] } split/;/,$class;
	my %transposonName;
	foreach my $l (@loc)
	{
		my $transposonStrand=($strand eq $$l[1])?"antisense":"sense";
		$transposonName{$$l[0]}=$transposonStrand;
	}
	my $c=1/scalar(keys %transposonName);
	push @{$item2coords{"$chrom;$strand;$_;$transposonName{$_}"}},[$start,$end,$sample,$c] foreach keys %transposonName; 
}
close IN;

my @results;
foreach my $item (keys %item2coords)
{
	my @sortedCoords=sort{ $a->[0]<=>$b->[0] } @{$item2coords{$item}};
	my ($chrom,$strand,$tName,$tStrand)=split(/;/,$item);
	my ($mergeStart,$mergeEnd,$mergeSample,$mergeCounts)=@{shift @sortedCoords};
	my %sampleCounts;
	$sampleCounts{$mergeSample}=$mergeCounts;
	foreach my $rangeRef (@sortedCoords) 
	{
    		my ($rangeStart,$rangeEnd,$rangeSample,$rangeCounts)=@{$rangeRef};
    		if($rangeEnd<=$mergeEnd)
		{
			$sampleCounts{$rangeSample}+=$rangeCounts;
			next;
		}
		if($rangeStart>=$mergeEnd)
		{
			my $count="";
			$count.=$_.",".$sampleCounts{$_}.";" foreach keys %sampleCounts;
			push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$count,$strand,$tStrand];
			($mergeStart,$mergeEnd,$mergeSample,$mergeCounts)=($rangeStart,$rangeEnd,$rangeSample,$rangeCounts);
			%sampleCounts=();
			$sampleCounts{$mergeSample}=$mergeCounts;
		}
		else
		{
			$mergeEnd=$rangeEnd;
			$sampleCounts{$rangeSample}+=$rangeCounts;
		}
	}
	my $count="";
	$count.=$_.",".$sampleCounts{$_}.";" foreach keys %sampleCounts;
	push @results,[$chrom,$mergeStart,$mergeEnd,$tName,$count,$strand,$tStrand] if $mergeEnd;
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
