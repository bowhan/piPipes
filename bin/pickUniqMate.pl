#!/share/bin/perl
use List::Util qw(sum);
use Bio::Seq;

die "perl $0 <mate sam with header> <uniq bed>\n" if @ARGV<1;

open in,$ARGV[1];
my %uniq;
while(<in>)
{
	chomp;
	my @f=split;
	$uniq{$f[3]}=[@f];
}
close in;

open in,$ARGV[0];
my (%te,@ref,%ref);
while(<in>)
{
	chomp;
	my @f=split/\t/,$_,12;
	# headers
	if(/^\@SQ/)
	{
		my ($sn,$ln)=/SN:(.*?)\tLN:(\d+)/;
		push @ref,[$sn,$ln];
		$ref{$sn}=$#ref;
		next;
	}

	# unmapped
	next if $f[2] eq "*";
	
	# alignments
	if($f[11]=~/XT:A:/)
	{
		my ($rnum)=$f[1]=~/(\d)$/;
		# CIGAR
		my (@cigar_m)=$f[5]=~/(\d+)M/g;
		my (@cigar_d)=$f[5]=~/(\d+)D/g;
		my (@cigar_s)=$f[5]=~/(\d+)S/g;
		my (@cigar_i)=$f[5]=~/(\d+)I/g;
		my $aln_ln=sum(@cigar_m,@cigar_d);
		
		my $strand="+";
        	if($f[1]=~/r/)
        	{
                	my $seq=Bio::Seq->new(-seq=>$f[9]);
                	$f[9]=$seq->revcom->seq;
                	$strand="-";
        	}

		# align to the junctions
		if(($f[3]+$aln_ln-1)>${$ref[$ref{$f[2]}]}[1])
		{
			if(($f[3]+($aln_ln-1)/2)>${$ref[$ref{$f[2]}]}[1])
			{
				$f[2]=${$ref[$ref{$f[2]}+1]}[0];
				$f[3]=1;
				$aln_ln=$aln_ln-(${$ref[$ref{$f[2]}]}[1]-$f[3]+1);
			}
			else
			{
				$aln_ln=${$ref[$ref{$f[2]}]}[1]-$f[3]+1;
			}
		}

		$pe{$f[0]}{$rnum}=$f[2].",".$strand."$f[3]".";";

		# XA tag
		if($f[11]=~/XA:Z:/)
		{
			my ($xa)=$f[11]=~/XA:Z:(.*);$/; 
			my @xa=split(";",$xa);
			$pe{$f[0]}{$rnum}.=join(",",(split/,/)[0,1]).";" foreach @xa;
		}
	}
}
close in;

foreach my $id (keys %pe)
{
	next if exists $pe{$id}{1} && exists $pe{$id}{2} && exists $uniq{$id."/1"} && exists $uniq{$id."/2"};
	foreach my $rid (keys %{$pe{$id}})
	{
		my $mate_id=($rid==1)?2:1;
		if(exists $uniq{$id."/".$mate_id})
		{
			${$uniq{$id."/".$mate_id}}[4]=$pe{$id}{$rid};
			print join("\t",@{$uniq{$id."/".$mate_id}}),"\n";
		}
	}
}
