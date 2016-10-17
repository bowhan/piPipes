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

# check whether the required enviromental variables have been exported correctly
declare -a ExportedVariableNames=(AnnotationDir GenomeAllmapBed2 PdfDir StatsTable CurrentInputNamePrefix GenomicHairpinReads)
for var in ${ExportedVariableNames[@]}; do
	if [[ -z ${!var} ]]; then echo2 "env variable $var is not exported, make sure you run the piPipes as a whole" error; fi
done
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
_summary_prefix=$1 # prefix to store the summary file
_smRNA_sum=${_summary_prefix}.smRNA.sum # summary for total small RNA
_siRNA_sum=${_summary_prefix}.siRNA.sum # summary for siRNA
_piRNA_sum=${_summary_prefix}.piRNA.sum # summary for piRNA
_intersect_outdir=$2 # output directory to store files
_seed=${RANDOM}${RANDOM}${RANDOM}${RANDOM} # random name
if [[ ! -s $AnnotationDir/genomic_features ]]; then 
	echo2 "Missing or empty $AnnotationDir/genomic_features file, cannot proceed with genomic structure intersecting..." error
fi
. $AnnotationDir/genomic_features # reading the information to intersect with, as well as some other annotation files
_all_bed=`basename ${GenomeAllmapBed2%bed2}x_rpmk_MASK.bed2` # names for the file genernated here

# get rid of tRNA, rRNA, snoRNA, if the genomic_features file defines MASK variable
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.intersect_remove_mask ]]; then
	if [[ -z $MASK ]]; then
		echo2 "undefined \$MASK. No masking will be done. If you would like to mask rRNA, tRNA, et al., please edit the $AnnotationDir/genomic_features file." "warning"
		if [[ ! -s $_intersect_outdir/${_all_bed} ]]; then ln -s $GenomeAllmapBed2 $_intersect_outdir/${_all_bed}; fi
		awk '{total[$7]=$4; if ($5==1) {unique_reads+=$4; ++unique_species}}END{for (seq in total) {all_reads+=total[seq]; ++all_species}; printf "%d\t%d\t%d\t%d\t", unique_reads, all_reads, unique_species, all_species}' \
			$_intersect_outdir/${_all_bed} \
			> $_intersect_outdir/.stats \
		&& touch $SentinelDir/${RunUid}.status.${Step}.intersect_remove_mask
	else
		bedtools_piPipes intersect -v -wa -a $GenomeAllmapBed2 -b $MASK \
		| tee $_intersect_outdir/${_all_bed} \
		| awk '{total[$7]=$4; if ($5==1) {unique_reads+=$4; ++unique_species}}END{for (seq in total) {all_reads+=total[seq]; ++all_species}; printf "%d\t%d\t%d\t%d\t", unique_reads, all_reads, unique_species, all_species}' \
		> $_intersect_outdir/.stats \
		&& touch $SentinelDir/${RunUid}.status.${Step}.intersect_remove_mask
	fi
else 
	echo2 "Removing masked feature with intersection has been done previously" warning
fi

print_header $_smRNA_sum
print_header $_siRNA_sum
print_header $_piRNA_sum

# doing intersecting and counting
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.intersect_all_features ]]; then
	para_file=$_intersect_outdir/${_seed}.intersect.para
	for t in ${TARGETS[@]}; do
		if [[ -z ${!t} ]]; then echo2 "variable ${t} was not defined, please exam your $AnnotationDir/genomic_features file" error; fi
		if [[ ! -f ${!t} ]]; then echo2 "${!t} cannot be found for ${t}, please exam your $AnnotationDir/genomic_features file" error; fi
		echo "bash piPipes_smallRNA_intersect.sh $_intersect_outdir/${_all_bed} ${t} ${!t} $_intersect_outdir/.stats" \
			>> $para_file
	done
	if [[ -s $para_file ]]; then
		ParaFly -c $para_file -CPU $Threads -failed_cmds ${para_file}.failedCommands &> /dev/null \
		&& rm -rf ${para_file}* \
		&& touch $SentinelDir/${RunUid}.status.${Step}.intersect_all_features \
		|| echo2 "Not all jobs in ${para_file} has completed" warning
	else 
		echo2 "No need to run intersecting for any features"
	fi
else 
	echo2 "Intersecting with features have been done previously" warning
fi

_pdfs=""
for t in ${TARGETS[@]}; do 
	if [[ -s $_intersect_outdir/${_all_bed}.intersect_with_${t}.unique_species.pdf ]]; then 
		_pdfs=${_pdfs}" "$_intersect_outdir/${_all_bed}.intersect_with_${t}.unique_species.pdf
	fi

	if [[ -s $_intersect_outdir/${_all_bed}.intersect_with_${t}.all_reads.pdf ]]; then 
		_pdfs=${_pdfs}" "$_intersect_outdir/${_all_bed}.intersect_with_${t}.all_reads.pdf
	fi 

	if [[ -s $_intersect_outdir/${_all_bed}.intersect_with_${t}.3p.pdf ]]; then 
		_pdfs=${_pdfs}" "$_intersect_outdir/${_all_bed}.intersect_with_${t}.3p.pdf
	fi 

	echo -ne "${t}\t" >> $_smRNA_sum
	cat $_intersect_outdir/.stats.${t}.smRNA >> $_smRNA_sum
	echo -ne "${t}\t" >> $_siRNA_sum
	cat $_intersect_outdir/.stats.${t}.siRNA >> $_siRNA_sum
	echo -ne "${t}\t" >> $_piRNA_sum
	cat $_intersect_outdir/.stats.${t}.piRNA >> $_piRNA_sum
done

# draw pie chart from exclusive groups
if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.pie ]]; then
	echo2 "Making pie plot to mutually exclusive gruops"
	rm -rf $_intersect_outdir/exclusive_genomic_feature.bed ${StatsTable}.exclusive_genomic_feature.count
	echo -ne "miRNAs\t" > ${StatsTable}.exclusive_genomic_feature.count 
	grep 'miRNA_hairpin_reads' $StatsTable \
		| awk '{print $NF}' >> ${StatsTable}.exclusive_genomic_feature.count \
	&& echo "miRNAs" > ${StatsTable}.exclusive_genomic_feature.order

	if [[ ! -z ${TARGETS_EXCLUSIVE[@]+"${TARGETS_EXCLUSIVE[@]}"} ]]; then
		for t in ${TARGETS_EXCLUSIVE[@]}; do
			# we run this everytime to ensure modification on the files
			echo "$t" >> ${StatsTable}.exclusive_genomic_feature.order
			bedSort ${!t} /dev/stdout \
			| bedtools_piPipes merge -i stdin \
			| awk -v name=$t 'BEGIN{OFS="\t"}{print $1,$2,$3,name}' \
				>> $_intersect_outdir/exclusive_genomic_feature.bed \
			|| echo2 "Failed to sort and merge exclude_genomic_feature" error
		done

		if [[ -f $_intersect_outdir/exclusive_genomic_feature.bed ]]; then
			bedtools_piPipes intersect -wo -a $GenomeAllmapBed2 -b $_intersect_outdir/exclusive_genomic_feature.bed \
				> ${StatsTable}.exclusive_genomic_feature.bed \
			&& awk '{ct[$(NF-2)]+=$4/$5/$NF}END{for (f in ct) {print f"\t"ct[f]}}' \
				${StatsTable}.exclusive_genomic_feature.bed \
				>> ${StatsTable}.exclusive_genomic_feature.count \
			&& awk -v total=$GenomicHairpinReads 'BEGIN{OFS="\t"}{ \
				if (ARGIND==1) {\
					a+=$2; \
					b[$1]=$2;\
				} else {\
					printf "%s\t%.1f\n",$1, (b[$1]?b[$1]:0);\
				} \
			} END{ \
				printf "unannotated\t%.1f\n", total-a \
			}' \
				${StatsTable}.exclusive_genomic_feature.count \
				${StatsTable}.exclusive_genomic_feature.order \
				> ${StatsTable}.exclusive_genomic_feature.count1 \
			&& mv ${StatsTable}.exclusive_genomic_feature.count1 ${StatsTable}.exclusive_genomic_feature.count \
			&& Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_pie.R \
				$PdfDir/${CurrentInputNamePrefix}.pie \
				${StatsTable}.exclusive_genomic_feature.count \
			&& awk -v siRNA_bot=$siRNA_bot -v siRNA_top=$siRNA_top '{l=$3-$2; if (l>=siRNA_bot && l<=siRNA_top) ct[$(NF-2)]+=$4/$5/$NF}END{for (f in ct) {print f"\t"ct[f]}}' \
				${StatsTable}.exclusive_genomic_feature.bed \
				> ${StatsTable}.exclusive_genomic_feature.siRNA.count \
			&& gawk 'BEGIN{OFS="\t"}{ if (ARGIND==1) {a+=$2; b[$1]=$2;} else {printf "%s\t%.1f\n",$1, (b[$1]?b[$1]:0);}}' \
				${StatsTable}.exclusive_genomic_feature.siRNA.count \
				${StatsTable}.exclusive_genomic_feature.order \
				> ${StatsTable}.exclusive_genomic_feature.siRNA.count1 \
			&& mv ${StatsTable}.exclusive_genomic_feature.siRNA.count1 ${StatsTable}.exclusive_genomic_feature.siRNA.count \
			&& Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_pie.R \
				$PdfDir/${CurrentInputNamePrefix}.siRNA.pie \
				${StatsTable}.exclusive_genomic_feature.siRNA.count \
			&& awk -v piRNA_bot=$piRNA_bot -v piRNA_top=$piRNA_top '{l=$3-$2; if (l>=piRNA_bot && l<=piRNA_top) ct[$(NF-2)]+=$4/$5/$NF}END{for (f in ct) {print f"\t"ct[f]}}' \
				${StatsTable}.exclusive_genomic_feature.bed \
				> ${StatsTable}.exclusive_genomic_feature.piRNA.count \
			&& gawk 'BEGIN{OFS="\t"}{ if (ARGIND==1) {a+=$2; b[$1]=$2;} else {printf "%s\t%.1f\n",$1, (b[$1]?b[$1]:0);}}' \
				${StatsTable}.exclusive_genomic_feature.piRNA.count \
				${StatsTable}.exclusive_genomic_feature.order \
				> ${StatsTable}.exclusive_genomic_feature.piRNA.count1 \
			&& mv ${StatsTable}.exclusive_genomic_feature.piRNA.count1 ${StatsTable}.exclusive_genomic_feature.piRNA.count \
			&& Rscript --slave ${PIPELINE_DIRECTORY}/bin/piPipes_draw_pie.R \
				$PdfDir/${CurrentInputNamePrefix}.piRNA.pie \
				${StatsTable}.exclusive_genomic_feature.piRNA.count \
			&& rm -rf $_intersect_outdir/exclusive_genomic_feature.bed ${StatsTable}.exclusive_genomic_feature.order ${StatsTable}.exclusive_genomic_feature.bed \
			&& touch $SentinelDir/${RunUid}.status.${Step}.pie \
			|| echo2 "failed to generate pie plot"
		fi 
	fi
else 
	echo2 "Pie plot has been made" warning
fi

if [[ ! -f $SentinelDir/${RunUid}.status.${Step}.combine_feature_pdf \
   && ! -z ${_pdfs} ]]; then
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$PdfDir/${CurrentInputNamePrefix}.features.pdf ${_pdfs} \
	&& rm -rf ${_pdfs} \
	&& touch $SentinelDir/${RunUid}.status.${Step}.combine_feature_pdf \
	|| echo2 "Failed to merge pdf from features intersecting... check gs... Or use your favorarite pdf merge tool by editing line$LINENO in $0" warning
fi