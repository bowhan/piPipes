# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) 2014  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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
GENOME_ALLMAP_BED2=`readlink -f $1` # genome wide alignment of all mappers, in bed2 format
SUMMARY_PREFIX=$2 # prefix to store the summary file
smRNA_SUM=${SUMMARY_PREFIX}.smRNA.sum # summary for total small RNA
siRNA_SUM=${SUMMARY_PREFIX}.siRNA.sum # summary for siRNA
piRNA_SUM=${SUMMARY_PREFIX}.piRNA.sum # summary for piRNA
CPU=$3 # CPU to use
INTERSECT_OUTDIR=$4 # output directory to store files
SEED=${RANDOM}${RANDOM}${RANDOM}${RANDOM} # random name
[ ! -s $COMMON_FOLDER/genomic_features ] && echo2 "Missing or empty $COMMON_FOLDER/genomic_features file, cannot proceed with genomic struture intersecting..." "error"
. $COMMON_FOLDER/genomic_features # reading the information to intersect with, as well as some other annotation files
ALL_BED=`basename ${GENOME_ALLMAP_BED2%bed2}x_rpmk_MASK.bed2` # names for the file genernated here
# get rid of tRNA, rRNA, snoRNA...
if [ -z $MASK ]; then
    echo2 "undefined \$MASK. No masking will be done. If you would like to mask rRNA, tRNA, et al., please edit the $COMMON_FOLDER/genomic_features file." "warning"
    [ ! -s $INTERSECT_OUTDIR/${ALL_BED} ] && ln -s $GENOME_ALLMAP_BED2 $INTERSECT_OUTDIR/${ALL_BED}
    awk '{total[$7]=$4; if ($5==1) {unique_reads+=$4; ++unique_species}}END{for (seq in total) {all_reads+=total[seq]; ++all_species}; printf "%d\t%d\t%d\t%d\t", unique_reads, all_reads, unique_species, all_species}' $INTERSECT_OUTDIR/${ALL_BED} > $INTERSECT_OUTDIR/.stats
else
    bedtools_piPipes intersect -v -wa -a $GENOME_ALLMAP_BED2 -b $MASK | \
    tee $INTERSECT_OUTDIR/${ALL_BED} | \
    awk '{total[$7]=$4; if ($5==1) {unique_reads+=$4; ++unique_species}}END{for (seq in total) {all_reads+=total[seq]; ++all_species}; printf "%d\t%d\t%d\t%d\t", unique_reads, all_reads, unique_species, all_species}' > $INTERSECT_OUTDIR/.stats
fi

print_header $smRNA_SUM
print_header $siRNA_SUM
print_header $piRNA_SUM

# doing intersecting and counting
para_file=$INTERSECT_OUTDIR/${SEED}.intersect.para
for t in ${TARGETS[@]}
do \
	echo "bash $DEBUG piPipes_smallRNA_intersect.sh $INTERSECT_OUTDIR/${ALL_BED}  ${t} ${!t} $INTERSECT_OUTDIR/.stats" >> $para_file
done
ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>&2 && \
rm -rf ${para_file}*

PDFs=""
for t in ${TARGETS[@]}
do \
	[ -s $INTERSECT_OUTDIR/${ALL_BED}.intersect_with_${t}.unique_species.pdf ] && PDFs=${PDFs}" "$INTERSECT_OUTDIR/${ALL_BED}.intersect_with_${t}.unique_species.pdf
	[ -s $INTERSECT_OUTDIR/${ALL_BED}.intersect_with_${t}.all_reads.pdf ] && PDFs=${PDFs}" "$INTERSECT_OUTDIR/${ALL_BED}.intersect_with_${t}.all_reads.pdf
	echo -ne "${t}\t" >> $smRNA_SUM
	cat $INTERSECT_OUTDIR/.stats.${t}.smRNA >> $smRNA_SUM
	echo -ne "${t}\t" >> $siRNA_SUM
	cat $INTERSECT_OUTDIR/.stats.${t}.siRNA >> $siRNA_SUM
	echo -ne "${t}\t" >> $piRNA_SUM
	cat $INTERSECT_OUTDIR/.stats.${t}.piRNA >> $piRNA_SUM
done

# draw pie chart from exclusive groups
rm -rf $INTERSECT_OUTDIR/exclusive_genomic_feature.bed
echo -ne "miRNAs\t" > ${TABLE}.exclusive_genomic_feature.count && \
grep 'miRNA hairpin reads' $TABLE | cut -f2 >> ${TABLE}.exclusive_genomic_feature.count && \
echo "miRNAs" > ${TABLE}.exclusive_genomic_feature.order
for t in ${TARGETS_EXCLUSIVE[@]}
do
	# we run this everytime to ensure modification on the files
	echo "$t" >> ${TABLE}.exclusive_genomic_feature.order
	bedSort ${!t} /dev/stdout | bedtools_piPipes merge -i stdin | awk -v name=$t 'BEGIN{OFS="\t"}{print $1,$2,$3,name}' >> $INTERSECT_OUTDIR/exclusive_genomic_feature.bed
done

[ -f $INTERSECT_OUTDIR/exclusive_genomic_feature.bed ] && \
	bedtools_piPipes intersect -wo -a $GENOME_ALLMAP_BED2 -b $INTERSECT_OUTDIR/exclusive_genomic_feature.bed | \
		awk '{ct[$(NF-2)]+=$4/$5/$NF}END{for (f in ct) {print f"\t"ct[f]}}' >> ${TABLE}.exclusive_genomic_feature.count && \
	awk -v total=$TOTAL_GENOME_MAPPING_READS 'BEGIN{OFS="\t"}{ if (ARGIND==1) {a+=$2; b[$1]=$2;} else {printf "%s\t%.1f\n",$1, (b[$1]?b[$1]:0);}}END{printf "unannotated\t%.1f\n", total-a }' ${TABLE}.exclusive_genomic_feature.count ${TABLE}.exclusive_genomic_feature.order > ${TABLE}.exclusive_genomic_feature.count1 && \
	mv ${TABLE}.exclusive_genomic_feature.count1 ${TABLE}.exclusive_genomic_feature.count && \
	Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_pie.R $PDF_DIR/${PREFIX}.pie ${TABLE}.exclusive_genomic_feature.count && \
	rm -rf $INTERSECT_OUTDIR/exclusive_genomic_feature.bed ${TABLE}.exclusive_genomic_feature.order

( gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PDF_DIR/${PREFIX}.features.pdf ${PDFs} && rm -rf ${PDFs} ) || \
echo2 "Failed to merge pdf from features intersecting... check gs... Or use your favorarite pdf merge tool by editing line$LINENO in $0" "warning"
