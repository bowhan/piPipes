#!/share/bin/perl
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <*.excision.cluster.rpmk.refined.bp>\n" if @ARGV<0;

my $title=$ARGV[0];
if ($title =~ /annotation/) {
    $title =~ s/excision.cluster.annotation.refined.bp/sorted.bam/;
}
else {$title =~ s/excision.cluster.rpmk.refined.bp/sorted.bam/;}
#system("samtools index /home/wangj2/scratch/bill/bill_genomic/$title");

my %chrs=();
system("samtools view -H $title > header");
open (input, "<header") or die "Can't open header since $!\n";
while (my $line=<input>) {
    if ($line =~ /^\@SQ/) {
	my @a=split(/\t/, $line);
	for my $j (0..$#a) {
            if ($a[$j] =~ /^SN:/) {
		$a[$j] =~ s/^SN://;
		$chrs{$a[$j]}=1;
            }
	}
    }
}
close input;
system("rm header");

open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
my $header=<input>;
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);

    my $left=0;
    my $right=0;
    if ($a[4] eq "") {$left=$a[1];}
    else {
	my @t=split(/\,/, $a[4]);
	my @p=split(/\(/, $t[$#t]);
	$left=$p[0];
    }
    if ($a[5] eq "") {$right=$a[2];}
    else {
	my @t=split(/\,/, $a[5]);
	my @p=split(/\(/, $t[0]);
	$right=$p[0];
    }

    my $leftlower=$left-500;
    my $leftupper=$left+500;
    my $rightlower=$right-500;
    my $rightupper=$right+500;
    my $chr_num=$a[0];
    $chr_num =~ s/chr//;
    if (($chrs{$a[0]} == 1) && (! defined $chrs{$chr_num})) {$chr_num=$a[0];}
    system("samtools view -bu $title $chr_num\:$leftlower\-$leftupper $chr_num\:$rightlower\-$rightupper > temp.bam");
    system("samtools view -Xf 0x2 temp.bam > temp.sam");
    
    open in,"temp.sam";
    my %ps=();
    my %me=();
    my %uniqp=();
    my %uniqm=();
    my $ref_sup=0;

    while(<in>)
    {
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;
	
	## XT:A:* 
	my ($xt)=$f[11]=~/XT:A:(.)/;
	
	## Coordinate
	my $coor=$f[3];
	if ($f[1]=~/r/)
	{
	    if ($xt eq "U") {$uniqm{$f[0]}=1;}
	    my (@cigar_m)=$f[5]=~/(\d+)M/g;
	    my (@cigar_d)=$f[5]=~/(\d+)D/g;
	    my (@cigar_s)=$f[5]=~/(\d+)S/g;
	    my (@cigar_i)=$f[5]=~/(\d+)I/g;
	    my $aln_ln=sum(@cigar_m,@cigar_d);
	    $me{$f[0]}=$f[3]+$aln_ln-1;
	}
	elsif ($f[1]=~/R/) {
	    $ps{$f[0]}=$f[3];
	    if ($xt eq "U") {$uniqp{$f[0]}=1;}
	}
	
#	${$pe{$f[0]}}[$rnum-1]=[$xt,$coor];
    }
    close in;

    foreach my $id (keys %ps)
    {
#	my @rid=@{$pe{$id}};
	
#	if(($rid[0][0] eq "U" && $rid[1][0] eq "M") || ($rid[0][0] eq "M" && $rid[1][0] eq "U"))
#	{
#	    $soft_clip++;
#	    print "$id\n";
#	}
	
	if ((defined $me{$id})&&((defined $uniqp{$id})||(defined $uniqm{$id})))
	{	    
            if (((($ps{$id}+5)<=$right)&&($me{$id}>$right)&&($uniqm{$id}==1)) || ((($me{$id}-5)>=$left)&&($ps{$id}<$left)&&($uniqp{$id}==1))) {	    
		$ref_sup++;
	    }
	}
    }

    print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$left\t$right\t$ref_sup\n";
    system("rm temp.sam temp.bam");
}

