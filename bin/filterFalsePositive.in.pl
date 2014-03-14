#!/share/bin/perl
use List::Util qw(max min);
#system("windowBed -a $ARGV[0] -b /home/wangj2/flycommon/all_transposons.dml.rmskCrossmatch.bed -sw -r 1000 -l 0 > tmp");
#system("windowBed -a $ARGV[0] -b /home/wangj2/flycommon/soo.trnalnpos.map2.sort.bed -sw -r 1000 -l 0 > tmp");
system("bedtools window -a $ARGV[0] -b $ARGV[1] -sw -r 1000 -l 0 > tmp");

open in,"tmp";
my %read;
while(<in>)
{
	chomp;
	split/\t/;

	## if the same tpye of transposons
	my @loc=map { [/(.*?),(\+|-)(.*)/] } split/;/,$_[4];
	foreach my $l (@loc)
	{
		if($$l[0] eq $_[9])
		{	
			## if the same strand of transposons
			if((($_[5] eq $$l[1]) && ($_[11] eq "-")) || (($_[5] ne $$l[1]) && ($_[11] eq "+")))
			{
				## if the fragments of the exists transposons
				{
					my $s=max($$l[2],$_[12]);
					my $e=min(($$l[2]+$_[2]-$_[1]),$_[13]);
					if($s<$e)
					{
						print join("\t",@_[0..5]),"\n" if not exists $read{$_[3]};
						$read{$_[3]}=1;
					}
				}
			}
		}
	}
}
close in;

system("rm tmp");
