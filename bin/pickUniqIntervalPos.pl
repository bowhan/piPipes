#!/share/bin/perl
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <sam> <fragment_size>\n" if @ARGV<1;
open in,$ARGV[0];
my %pe;
while(<in>)
{
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;

	## XT:A:* 
	my ($xt)=$f[11]=~/XT:A:(.)/;

	my $strand="+";

	## parse CIGAR
	if(($f[1]=~/R/)&&($f[8] > $ARGV[1])&&($f[8] <= 10000))
        {
                # CIGAR
                my (@cigar_m)=$f[5]=~/(\d+)M/g;
                my (@cigar_d)=$f[5]=~/(\d+)D/g;
                my (@cigar_s)=$f[5]=~/(\d+)S/g;
                my (@cigar_i)=$f[5]=~/(\d+)I/g;
                my $aln_ln=sum(@cigar_m,@cigar_d);
		
#		print $f[2],"\t",$f[3]-1+$aln_ln,"\t",$f[3]+$f[8],"\t$f[0]/$rnum\t","\n";
		if ($f[2] =~ /^\d{1,2}$/) {$f[2]="chr$f[2]";}
		print $f[2],"\t",$f[3]-6+$aln_ln,"\t",$f[7]+5,"\t$f[0]/$rnum\t","\n";
	}
}
close in;

