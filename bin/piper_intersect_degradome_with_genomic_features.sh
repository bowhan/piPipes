
# Degradome pipeline from piper: https://github.com/bowhan/piper.git
# piper: https://github.com/bowhan/piper.git
# An integrated pipeline for piRNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Weng lab
# University of Massachusetts Medical School

#####################
#function declartion#
#####################
function print_header {
	echo -ne "feature\t" > $1
	echo -ne "total_lib_all_mapper_reads\ttotal_feature_all_mapper_reads\tfeature_all_mapper_percentage\tfeature_sense_all_mapper_reads\tfeature_antisense_all_mapper_reads\tfeature_all_mapper_sense_fraction\t" >> $1
	echo -ne "total_lib_unique_mapper_reads\ttotal_feature_unique_mapper_reads\tfeature_unique_mapper_percentage\tfeature_sense_unique_mapper_reads\tfeature_antisense_unique_mapper_reads\tfeature_unique_mapper_sense_fraction\t" >> $1
	echo -ne "total_lib_unique_mapper_species\ttotal_feature_unique_mapper_species\tfeature_unique_mapper_percentage\tfeature_sense_unique_mapper_species\tfeature_antisense_unique_mapper_species\tfeature_unique_mapper_sense_fraction\n" >> $1
}
###############
#configuration#
###############
GENOME_ALLMAP_BED12=$1 # genome wide alignment of all mappers, in bed12 format, with $5 being 1/NTM
SUMMARY=$2 # summary file
CPU=$3 # CPU to use
INTERSECT_OUTDIR=$4 # output directory to store files
SEED=${RANDOM}${RANDOM}${RANDOM}${RANDOM} # random name
. $COMMON_FOLDER/genomic_features # reading the information to intersect with, as well as some other annotation files
ALL_BED=`basename ${GENOME_ALLMAP_BED12%bed*}x_rpmk_rtRNA.bed12` # names for the file genernated here
# get rid of tRNA, rRNA, snoRNA...
bedtools_piper intersect -v -wa -split -a $GENOME_ALLMAP_BED12 -b $rtRNA | tee $INTERSECT_OUTDIR/${ALL_BED} | awk '{all_reads+=$5; if ($5==1) {unique_reads+=$4; ++unique_species}}END{printf "%d\t%d\t%d\n", unique_reads, all_reads, unique_species;}' > $INTERSECT_OUTDIR/.stats
print_header $SUMMARY
# doing intersecting and counting
para_file=$INTERSECT_OUTDIR/${SEED}.intersect.para
for t in ${TARGETS[@]}
do \
	echo "bash $DEBUG piper_degradome_intersect.sh $INTERSECT_OUTDIR/${ALL_BED} ${t} ${!t} $INTERSECT_OUTDIR/.stats $SRA_UNIQ_BED2" >> $para_file
done
ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
rm -rf ${para_file}*

PDFs=""
for t in ${TARGETS[@]}
do \
	PDFs=${PDFs}" "$INTERSECT_OUTDIR/${ALL_BED}.intersect_with_${t}.pdf
	echo -ne "${t}\t" >> $SUMMARY
	cat $INTERSECT_OUTDIR/.stats.${t} >> $SUMMARY
done

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.features.pdf ${PDFs} && \
rm -rf ${PDFs} || \
echo2 "Failed to merge pdf from features intersecting... check gs... Or use your favorarite pdf merge tool by editing line$LINENO in $0" "warning"


