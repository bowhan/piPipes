
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
bed=$1
feature_name=$2
feature_bed=$3
total_lib_stats_file=$4
unique_reads=`cut -f1 $total_lib_stats_file`
all_reads=`cut -f2 $total_lib_stats_file`
unique_species=`cut -f3 $total_lib_stats_file`
export ext_len=30 # extend 30nt up/downstream

# intersect with the bed without strand information
bedtools_piPipes intersect -split -wo -f 0.5 -a $bed -b $feature_bed > ${bed}.intersect_with_${feature_name} && \
bedtools_piPipes intersect -split -wo -f 0.5 -a ${bed}.intersect_with_${feature_name} -b $feature_bed -s > ${bed}.intersect_with_${feature_name}.S && \
bedtools_piPipes intersect -split -wo -f 0.5 -a ${bed}.intersect_with_${feature_name} -b $feature_bed -S > ${bed}.intersect_with_${feature_name}.AS && \
awk -v lib_uniq_reads=$unique_reads \
	-v lib_all_reads=$all_reads \
	-v lib_unique_species=$unique_species \
	'BEGIN {FS=OFS="\t";} \
	{ \
		c=$5/$NF; \
		if ($6==$18) \
		{ \
			sense_all_counter+=c; \
			if ($5==int($5)) \
			{ \
				sense_uniq_counter+=c; \
				sense_uniq_species++; \
			}\
		} \
		else \
		{ \
			antisense_all_counter+=c; \
			if ($5==int($5)) \
			{ \
				antisense_uniq_counter+=c; \
				antisense_uniq_species++; \
			} \
		} \
	} \
	END\
	{\
		printf "%.0f\t%.0f\t%.3f\t%.0f\t%.0f\t%.3f\t%.0f\t%.0f\t%.3f\t%.0f\t%.0f\t%.3f\t%.0f\t%.0f\t%.3f\t%.0f\t%.0f\t%.3f\n",\
		lib_all_reads,\
		sense_all_counter+antisense_all_counter,\
		lib_all_reads!=0?((sense_all_counter+antisense_all_counter)/lib_all_reads):0,\
		sense_all_counter,\
		antisense_all_counter,\
		(sense_all_counter+antisense_all_counter)!=0? (sense_all_counter/(sense_all_counter+antisense_all_counter)):0,\
		lib_uniq_reads,\
		sense_uniq_counter+antisense_uniq_counter,\
		lib_uniq_reads!=0?((sense_uniq_counter+antisense_uniq_counter)/lib_uniq_reads):0,\
		sense_uniq_counter,\
		antisense_uniq_counter,\
		(sense_uniq_counter+antisense_uniq_counter)!=0?(sense_uniq_counter/(sense_uniq_counter+antisense_uniq_counter)):0,\
		lib_unique_species,\
		sense_uniq_species+antisense_uniq_species,\
		lib_unique_species!=0? ((sense_uniq_species+antisense_uniq_species)/lib_unique_species):0,\
		sense_uniq_species,\
		antisense_uniq_species,\
		(sense_uniq_species+antisense_uniq_species)!=0? (sense_uniq_species/(sense_uniq_species+antisense_uniq_species)):0;\
}' ${bed}.intersect_with_${feature_name} > ${total_lib_stats_file}.${feature_name}

awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$9])) {printed[$9]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}    | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.species.5end_60.percentage
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$9])) {printed[$9]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.S  | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.species.5end_60.percentage.S
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$9])) {printed[$9]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.AS | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.species.5end_60.percentage.AS

# if small RNA 
if [ -n "$5" ] ; then 
	bedtools_piPipes intersect        -wa -u -f 0.5 -a $5   -b $feature_bed > ${5}.intersect_with_${feature_name} && \
	bedtools_piPipes intersect -split -wa -u -f 0.5 -a $bed -b $feature_bed | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,int($5),1,$6}' > ${bed}.intersect_with_${feature_name} && \
	piPipes_local_ping_pong -a ${5}.intersect_with_${feature_name} -b ${bed}.intersect_with_${feature_name} -p 1 > ${bed}.intersect_with_${feature_name}.pp.`basename $5`
	Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_degradome_features.R \
		${bed}.intersect_with_${feature_name} \
		${ext_len} \
		${bed}.intersect_with_${feature_name}.species.5end_60.percentage \
		${bed}.intersect_with_${feature_name}.pp.`basename $5` \
		1>&2 && \
	rm -rf ${5}.intersect_with_${feature_name} ${bed}.intersect_with_${feature_name} ${bed}.intersect_with_${feature_name}.S ${bed}.intersect_with_${feature_name}.AS
else
	Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_degradome_features.R \
		${bed}.intersect_with_${feature_name} \
		${ext_len} \
		${bed}.intersect_with_${feature_name}.species.5end_60.percentage \
		1>&2 && \
	rm -rf ${bed}.intersect_with_${feature_name} ${bed}.intersect_with_${feature_name}.S ${bed}.intersect_with_${feature_name}.AS
fi

