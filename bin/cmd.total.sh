rm tmp
#for i in flamBGFM flamKGFM flamEmbryo flamHets flamTranshets
#for i in harwich wXh24 wXh14 wXh21 whXw14 whXw21
#for i in armiTranshetsOvary armiTranshetsSoma armiHetsOvary armiHetsSoma rhinoTranshetsOvary rhinoTranshetsSoma rhinoHetsOvary rhinoHetsSoma qinTranshetsOvary qinHetsOvary w1118Ovary w1118Soma orerOvary orerSoma orerEmbryo
#for i in w1118Ovary w1118Soma qinTranshetsOvary qinHetsOvary qintestTranshetsOvary qintestHetsOvary 
for i in harwich.ovary harwichG20.ovary W1.ovary W1G20.ovary wXh1g.ovary whXh3g.ovary whXh5g.ovary whXh7g.ovary whXw3g.ovary whXw5g.ovary whXw7g.ovary
#for i in W1.ovary wXh1g.ovary whXh3g.ovary whXh5g.ovary whXw3g.ovary whXw5g.ovary armiHets.ovary armiHets.carcass armiTranshets.ovary armiTranshets.carcass flamBGFM.ovary flam.embryo flamHets.ovary flamKGFM.ovary flamTranshets.ovary harwich.ovary introgression2.ovary introgression2X3.ovary introgression3.ovary orer.embryo orer.ovary orer.carcass qinDf.ovary qinHets.ovary qinTMB.ovary qinTranshets.ovary rhinoHets.ovary rhinoHets.carcass rhinoTranshets.ovary rhinoTranshets.carcass w1.ovary w1.carcass whXw14d.ovary whXw21d.ovary wXh14d.ovary wXh21d.ovary wXh2_4d.ovary

do
	
	awk -F "\t" -v sample=$i '{OFS="\t"; print $1,$2,$3,sample,$5,$6}' $i.downsample.bam.unpair.uniq.transposons.filtered.bed >> tmp
	
	## Filter BS A{36}
	grep FBgn0000224_BS tmp | egrep "\+51|\-51" > tmp.BS

	## Merge Stalker
	ediff tmp diff tmp.BS > tmp2 	

done

perl /home/wangj2/jpp_findTransposonJumping/mergeTagsWithoutGap.pl tmp2 > dysgenic.uniq.transposons.filtered.woGap.bed
perl /home/wangj2/jpp_findTransposonJumping/mergeTagsWithGap.pl dysgenic.uniq.transposons.filtered.woGap.bed 500 > dysgenic.uniq.transposons.filtered.wGap.bed

rm tmp2 tmp.BS tmp

perl get_class.pl dysgenic.uniq.transposons.filtered.wGap.bed > dysgenic.uniq.transposons.filtered.wGap.class.bed
