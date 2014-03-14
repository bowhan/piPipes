#!/share/bin/perl
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <*.excision.cluster.rpmk> <Reference.2bit>\n" if @ARGV<1;

my $title=$ARGV[0];
if ($title =~ /annotation/) {
    $title =~ s/excision.cluster.annotation/sorted.bam/;
}
else {$title =~ s/excision.cluster.rpmk/sorted.bam/;}

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
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\s+/, $line);

    my $lower=$a[3]-100;
    my $upper=$a[4]+100;
    my $chr_num=$a[2];
    $chr_num =~ s/chr//;
    if (($chrs{$a[2]} == 1) && (! defined $chrs{$chr_num})) {$chr_num=$a[2];}
    system("samtools view -bu $title $chr_num\:$lower\-$upper > temp.bam");
    system("samtools view -Xf 0x2 temp.bam > temp.sam");

    my $leftseq="";
    my $rightseq="";

    my $ll=$a[3]-150;
    my $lu=$a[3]+150;
    system("twoBitToFa $ARGV[1] -seq=$a[2] -start=$ll -end=$lu left.fa");
    open (seq, "<left.fa") or die "Can't open left.fa since $!\n";
    my $head=<seq>;
    for my $k (0..5) {
	$head=<seq>;
	chomp($head);
	$leftseq=$leftseq."$head";
    }
    $leftseq=uc($leftseq);
    close seq;
    system("rm left.fa");

    my $rl=$a[4]-150;
    my $ru=$a[4]+150;
    system("twoBitToFa $ARGV[1] -seq=$a[2] -start=$rl -end=$ru right.fa");
    open (seq, "<right.fa") or die "Can't open right.fa since $!\n";
    my $head=<seq>;
    for my $k (0..5) {
	$head=<seq>;
	chomp($head);
	$rightseq=$rightseq."$head";
    }
    $rightseq=uc($rightseq);
    close seq;
    system("rm right.fa");
    

    open in,"temp.sam";
    my %pe=();
    while(<in>)
    {
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;
	
	## XT:A:* 
	my ($xt)=$f[11]=~/XT:A:(.)/;

	my $CIGAR=$f[5];
	$CIGAR =~ s/S//g;
	if ($f[5]=~/S/) {
	
	    ## Coordinate
            my $coor=-10;
	    my $overcoor=-10;
            my $strand="";
	    my @z=split(/M/, $f[5]);

            if (($f[5]=~/S$/)&&($f[1]=~/r/))
            {
		my (@cigar_m)=$f[5]=~/(\d+)M/g;
                my (@cigar_d)=$f[5]=~/(\d+)D/g;
                my (@cigar_s)=$f[5]=~/(\d+)S/g;
                my (@cigar_i)=$f[5]=~/(\d+)I/g;
                my $aln_ln=sum(@cigar_m,@cigar_d);
		$coor=$f[3]+$aln_ln-1;
                $strand="-";

		my (@clipped)=$z[1]=~/(\d+)S/g;
		my $cliplen=sum(@clipped);
#		print "$f[0]\n";
#		print "$cliplen\t";
		if ($cliplen >= 10) {
		    my $clipseq=substr($f[9], length($f[9])-$cliplen, $cliplen);
		    $overcoor = index($rightseq, $clipseq);
#		    print "$clipseq\t$rightseq\t$overcoor\t";
		    if ($overcoor > -1) {$overcoor += ($a[4] - 149);}
		}
#		print "\n";
            }
            elsif (($f[1]=~/R/)&&($z[0]=~/S/)) {
		$coor=$f[3]; $strand="+";
		my (@clipped)=$z[0]=~/(\d+)S/g;
                my $cliplen=sum(@clipped);
#		print "$f[0]\n";
#		print "$cliplen\t";
		if ($cliplen >= 10) {
		    my $clipseq=substr($f[9], 0, $cliplen);
		    $overcoor = index($leftseq, $clipseq);
#		    print "$clipseq\t$leftseq\t$overcoor\t";
		    if ($overcoor > -1) {$overcoor += ($a[3] - 150 + $cliplen);}
		}
#		print "\n";
	    }

            if ($coor > 0) {
                my $final="";
		if ($overcoor > 0) {
		    if ($strand eq "-") {$final="$coor\-$overcoor"."\($strand\)";}
		    else {$final="$overcoor\-$coor"."\($strand\)";}
		}
		else {$final=$coor."\($strand\)";}
		if (defined $pe{$final}) {$pe{$final}++;}
		else {$pe{$final}=1;}
            }

	}
    }
    close in;

    my $clip_site="";
    
    foreach my $coor (keys %pe)
    {
	if ($pe{$coor} >= 2) {
	    $clip_site=$clip_site."$coor\:$pe{$coor}\;";
	}
    }

    chop($clip_site);
    print "$a[2]\t$a[3]\t$a[4]\t$a[5]\t$clip_site\n";
    system("rm temp.sam temp.bam");
#    last;
}

