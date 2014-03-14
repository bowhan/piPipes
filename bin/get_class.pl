#!/share/bin/perl

# chr2L   19384   20049   FBgn0001283_jockey      wXh24,-1;harwich,-8,+5;whXw21,-4,+1;whXw14,-5;wXh14,+2; +-      sense   chr2L:19562.19645

my @sample=("$ARGV[1]");

print "chr\tstart\tend\ttransposonName\tstrand\ttransposonStrand\tbreak\tclass";
print "\t$_\_class\t$_\_plus\t$_\_minus" foreach @sample;
print "\n";
open in,$ARGV[0];
while(<in>)
{
	chomp;
	my($chrom,$start,$end,$transposonName,$class,$strand,$transposonStrand,$break)=split/\t/;
	my %classCounts;
	my ($tcplus,$tcminus)=(0,0);
	foreach $s (split/;/,$class)
	{
		my ($name,@counts)=split/,/,$s;
		foreach my $c (@counts)
		{
			my $strand=($c>0)?"+":"-";
			$classCounts{$name}{$strand}=$c;
			$tcplus+=$c if $c>0;
			$tcminus+=$c if $c<0;
		}
	}
	print "$chrom\t$start\t$end\t$transposonName\t$strand\t$transposonStrand\t$break";
	print "\t1p1" if $tcplus>0 && $tcminus<0;
	print "\t2p" if ($tcplus>1 && $tcminus==0) || ($tcplus==0 && $tcminus<-1);
	print "\tsingleton" if ($tcplus<=1 && $tcminus==0 && $tcplus>0) || ($tcplus==0 && $tcminus>=-1 && $tcminus<0);
	print "\tNone" if ($tcminus==0 && $tcplus==0);
	foreach my $s (@sample)
	{
		$classCounts{$s}{"+"}=0 if not exists $classCounts{$s}{"+"};
		$classCounts{$s}{"-"}=0 if not exists $classCounts{$s}{"-"};
		print "\t1p1" if $classCounts{$s}{"+"}>0 && $classCounts{$s}{"-"}<0;
		print "\t2p" if ($classCounts{$s}{"+"}>1 && $classCounts{$s}{"-"}==0) || ($classCounts{$s}{"+"}==0 && $classCounts{$s}{"-"}<-1);
		print "\tsingleton" if ($classCounts{$s}{"+"}<=1 && $classCounts{$s}{"-"}==0 && $classCounts{$s}{"+"}>0) || ($classCounts{$s}{"+"}==0 && $classCounts{$s}{"-"}>=-1 && $classCounts{$s}{"-"}<0);
		print "\tNone" if $classCounts{$s}{"+"}==0 && $classCounts{$s}{"-"}==0;
		print "\t",$classCounts{$s}{"+"},"\t",$classCounts{$s}{"-"};
	}
	print "\n";
}
close in;
