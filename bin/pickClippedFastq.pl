#!/share/bin/perl
use List::Util qw(sum);
use Bio::Seq;

die "perl $0 <input_prefix> <TE sequence database>\n" if @ARGV<1;

my %transposon_seq=();
my %transposon_revcom_seq=();
my $curr_seq="";
my $curr_transposon="";
open (input, "<$ARGV[1]") or die "Can't open $ARGV[1] since $!\n";
while (my $line=<input>) {
    chomp($line);
    if ($line =~ /^\>/) {
	if ($curr_transposon ne "") {
	    $transposon_seq{$curr_transposon}=uc($curr_seq);
	    my $seq=Bio::Seq->new(-seq=>$curr_seq);
	    $curr_seq=$seq->revcom->seq;
	    $transposon_revcom_seq{$curr_transposon}=uc($curr_seq);
	}
	my @a=split(/\s+/, $line);
	$a[0] =~ s/\>//;
	$curr_transposon=$a[0];
	$curr_seq="";
    }
    else {$curr_seq=$curr_seq.$line;}
}
close input;

open m1,">>$ARGV[0].clipped.reads.aln";

open (input, "<$ARGV[0].insertion.bp.bed") or die "Can't open $ARGV[0].insertion.bp.bed since $!\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);

    my $lower=$a[1]-15;
    my $upper=$a[2]+15;
    if (($lower > 0)&&($upper > 0))
    {
	system("samtools view -hX $ARGV[0].pair.bam $a[0]\:$lower\-$upper > temp.sam");
	
	open in,"temp.sam";
	my %pe1;
	my %pe2;
	while(<in>)
	{
	    chomp;
	    my @f=split/\t/,$_,12;
	    ## read number 1 or 2
	    my ($rnum)=$f[1]=~/(\d)$/;
	    
	    ## XT:A:* 
	    my ($xt)=$f[11]=~/XT:A:(.)/;
	    
	    if ($f[5]=~/S/) {
		
		## Coordinate                                                                                                                                    
		my $coor=-10;
		my $strand="";
		my $final="";
		my $clipseq="";
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
		    if ($cliplen >= 7) {
			$clipseq=substr($f[9], length($f[9])-$cliplen, $cliplen);
		    }
		}

                elsif (($f[1]=~/R/)&&($z[0]=~/S/))
                {
                    $coor=$f[3]; $strand="+";

                    my (@clipped)=$z[0]=~/(\d+)S/g;
                    my $cliplen=sum(@clipped);
                    if ($cliplen >= 7) {
                        $clipseq=substr($f[9], 0, $cliplen);
		    }
		}

		if ($clipseq ne "") {
		    my $flag=0;
		    while ((my $key, my $value) = each (%transposon_seq)) {
			my $seq=$value;
			if ($a[4] eq "antisense") {
			    $seq=$transposon_revcom_seq{$key};
			}
			if (($seq =~ /$clipseq/)&&($a[3] eq $key)) {
#			    print "$clipseq\n";
			    $final=$coor."\($strand\)";
			    if (defined $pe1{$final}) {
				if (length($clipseq) > length($pe1{$final})) {
				    $pe1{$final}=$clipseq;
				}
			    }
			    else {$pe1{$final}=$clipseq; $pe2{$final}=0;}
			    $flag=1;
			    last;
			}
		    }		    
		}
		
	    }	    
	}#while;
	close in;

        open in,"temp.sam";
        while(<in>)
        {
            chomp;
            my @f=split/\t/,$_,12;
            my ($rnum)=$f[1]=~/(\d)$/;
            my ($xt)=$f[11]=~/XT:A:(.)/;

            if ($f[5]=~/S/) {

                my $coor=-10;
                my $strand="";
		my $final="";
                my $clipseq="";
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
                    if ($cliplen >= 4) {
                        $clipseq=substr($f[9], length($f[9])-$cliplen, $cliplen);
                    }
                }

                elsif (($f[1]=~/R/)&&($z[0]=~/S/))
                {
                    $coor=$f[3]; $strand="+";

                    my (@clipped)=$z[0]=~/(\d+)S/g;
                    my $cliplen=sum(@clipped);
                    if ($cliplen >= 4) {
                        $clipseq=substr($f[9], 0, $cliplen);
                    }
                }

                if ($clipseq ne "") {
		    foreach my $coor (keys %pe1) {
			if (($coor =~ /\+/) && (substr($pe1{$coor}, length($pe1{$coor})-length($clipseq), length($clipseq)) eq $clipseq)) {
			    $pe2{$coor}++;
			}
			elsif (($coor =~ /\-/) && (substr($pe1{$coor}, 0, length($clipseq)) eq $clipseq)) {
			    $pe2{$coor}++;
			}
		    }
		}
	    }
	}
	close in;

	my $clip_site="";
	
	foreach my $coor (keys %pe2)
	{
	    $clip_site=$clip_site."$coor\:$pe2{$coor}\;";
	}
	chop($clip_site);
	print m1 "$a[0]\t$lower\t$upper\t$a[3]\t$clip_site\n";
	system("rm temp.sam");
    }
    else 
    {
	print m1 "$line\n";
    }
}
close input;
close m1;
