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
bed=$1 # geomic mapping in bed2 format
feature_name=$2 # name of this feature
feature_bed=$3 # genomic locus of this feature
total_lib_stats_file=$4 # statistics to write in this feature, with some information written (below)
unique_reads=`cut -f1 $total_lib_stats_file`
all_reads=`cut -f2 $total_lib_stats_file`
unique_species=`cut -f3 $total_lib_stats_file`
all_species=`cut -f4 $total_lib_stats_file`
export ext_len=30 # extend 30nt up/downstream // static const ...

# intersect with the bed without strand information, using a modified version of bedtools
bedtools_piPipes intersect -wo -f 0.5 -a $bed -b $feature_bed > ${bed}.intersect_with_${feature_name} && \
awk -v lib_uniq_reads=$unique_reads \
	-v lib_all_reads=$all_reads \
	-v lib_unique_species=$unique_species \
	-v lib_all_species=$all_species \
	'BEGIN {FS=OFS="\t";} \
	{ \
		c=$4/$5/$NF; \
		if ($6==$13) \
		{ \
			sense_all_counter+=c; \
			if ($5==1) \
			{ \
				sense_uniq_counter+=c; \
				sense_uniq_species++; \
			}\
		} \
		else \
		{ \
			antisense_all_counter+=c; \
			if ($5==1) \
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
}' ${bed}.intersect_with_${feature_name} > ${total_lib_stats_file}.${feature_name}.smRNA && \
awk -v l=$siRNA_bot -v h=$siRNA_top '$3-$2>=l&&$3-$2<=h' ${bed}.intersect_with_${feature_name} | tee ${bed}.intersect_with_${feature_name}.siRNA |  \
awk -v lib_uniq_reads=$unique_reads \
	-v lib_all_reads=$all_reads \
	-v lib_unique_species=$unique_species \
	-v lib_all_species=$all_species \
	'BEGIN {FS=OFS="\t";} \
	{ \
		c=$4/$5/$NF; \
		if ($6==$13) \
		{ \
			sense_all_counter+=c; \
			if ($5==1) \
			{ \
				sense_uniq_counter+=c; \
				sense_uniq_species++; \
			}\
		} \
		else \
		{ \
			antisense_all_counter+=c; \
			if ($5==1) \
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
}' > ${total_lib_stats_file}.${feature_name}.siRNA && \
awk -v l=$piRNA_bot -v h=$piRNA_top '$3-$2>=l&&$3-$2<=h' ${bed}.intersect_with_${feature_name} | tee ${bed}.intersect_with_${feature_name}.piRNA | \
awk -v lib_uniq_reads=$unique_reads \
	-v lib_all_reads=$all_reads \
	-v lib_unique_species=$unique_species \
	-v lib_all_species=$all_species \
	'BEGIN {FS=OFS="\t";} \
	{ \
		c=$4/$5/$NF; \
		if ($6==$13) \
		{ \
			sense_all_counter+=c; \
			if ($5==1) \
			{ \
				sense_uniq_counter+=c; \
				sense_uniq_species++; \
			}\
		} \
		else \
		{ \
			antisense_all_counter+=c; \
			if ($5==1) \
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
}' > ${total_lib_stats_file}.${feature_name}.piRNA && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	if ($5==1) \
	{ \
		l=$3-$2; \
		if (l>m) m=l; \
		if ( $6==$13 ) s[l]+=$4/$NF;\
		else as[l]+=$4/$NF; \
	} \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name} > ${bed}.intersect_with_${feature_name}.lendis && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	l=$3-$2; \
	if (l>m) m=l; \
	if ( $6==$13 ) s[l]+=$4/$5/$NF;\
	else as[l]+=$4/$5/$NF; \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name} > ${bed}.intersect_with_${feature_name}.allMapper.lendis && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	if ($5==1) \
	{ \
		l=$3-$2; \
		if (l>m) m=l; \
		if ( $6==$13 ) s[l]+=$4/$NF;\
		else as[l]+=$4/$NF; \
	} \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name}.siRNA > ${bed}.intersect_with_${feature_name}.siRNA.lendis && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	l=$3-$2; \
	if (l>m) m=l; \
	if ( $6==$13 ) s[l]+=$4/$5/$NF;\
	else as[l]+=$4/$5/$NF; \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name}.siRNA > ${bed}.intersect_with_${feature_name}.allMapper.siRNA.lendis && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	if ($5==1) \
	{ \
		l=$3-$2; \
		if (l>m) m=l; \
		if ( $6==$13 ) s[l]+=$4/$NF;\
		else as[l]+=$4/$NF; \
	} \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name}.piRNA > ${bed}.intersect_with_${feature_name}.piRNA.lendis && \
awk 'BEGIN{FS=OFS="\t"}\
{ \
	l=$3-$2; \
	if (l>m) m=l; \
	if ( $6==$13 ) s[l]+=$4/$5/$NF;\
	else as[l]+=$4/$5/$NF; \
}END\
{\
	for (d=1;d<=m;++d) \
	{\
		printf "%d\t%.0f\t%.0f\n", d, (s[d]?s[d]:0), (as[d]?as[d]:0); \
	}\
}' ${bed}.intersect_with_${feature_name}.piRNA > ${bed}.intersect_with_${feature_name}.allMapper.piRNA.lendis

# species
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name} | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.species.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.siRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.siRNA.species.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.piRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.piRNA.species.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name} | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.species.3end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.siRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.siRNA.species.3end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=1;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.piRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.piRNA.species.3end_60.percentage
# reads
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name} | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.reads.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.siRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.siRNA.reads.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="+") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.piRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.piRNA.reads.5end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name} | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.reads.3end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.siRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.siRNA.reads.3end_60.percentage && \
awk -v ext_len=$ext_len 'BEGIN{OFS="\t"} { if (($5==1)&&(!printed[$7])) {printed[$7]=1; if ($2>=ext_len) { for (i=1;i<=$4;++i) { if ($6=="-") { print $1,$2-ext_len,$2+ext_len+1,$4,$5,$6 } else { print $1,$3-ext_len-1,$3+ext_len,$4,$5,$6 }}}}}' ${bed}.intersect_with_${feature_name}.piRNA | bedtools_piPipes getfasta -fi $GENOME_FA -bed stdin -fo stdout -s -name -tab | piPipes_nuc_percentage.py $ext_len > ${bed}.intersect_with_${feature_name}.piRNA.reads.3end_60.percentage

# running ping-pong
bedtools_piPipes intersect -u -wa -f 0.5 -a $bed -b $feature_bed | tee ${bed}.intersect_with_${feature_name} | awk '$5==1' > ${bed}.intersect_with_${feature_name}.uniqueMappers
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name} -b ${bed}.intersect_with_${feature_name} -p 1 > ${bed}.intersect_with_${feature_name}.pp
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name}.uniqueMappers -b ${bed}.intersect_with_${feature_name}.uniqueMappers -p 1 > ${bed}.intersect_with_${feature_name}.uniqueMappers.pp

awk -v l=$siRNA_bot -v h=$siRNA_top '$3-$2>=l&&$3-$2<=h' ${bed}.intersect_with_${feature_name} | tee ${bed}.intersect_with_${feature_name}.siRNA | awk '$5==1' > ${bed}.intersect_with_${feature_name}.siRNA.uniqueMappers
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name}.siRNA -b ${bed}.intersect_with_${feature_name}.siRNA -p 1 > ${bed}.intersect_with_${feature_name}.siRNA.pp
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name}.siRNA.uniqueMappers -b ${bed}.intersect_with_${feature_name}.siRNA.uniqueMappers -p 1 > ${bed}.intersect_with_${feature_name}.siRNA.uniqueMappers.pp

awk -v l=$piRNA_bot -v h=$piRNA_top '$3-$2>=l&&$3-$2<=h' ${bed}.intersect_with_${feature_name} | tee ${bed}.intersect_with_${feature_name}.piRNA | awk '$5==1' > ${bed}.intersect_with_${feature_name}.piRNA.uniqueMappers
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name}.piRNA -b ${bed}.intersect_with_${feature_name}.piRNA -p 1 > ${bed}.intersect_with_${feature_name}.piRNA.pp
piPipes_local_ping_pong -a ${bed}.intersect_with_${feature_name}.piRNA.uniqueMappers -b ${bed}.intersect_with_${feature_name}.piRNA.uniqueMappers -p 1 > ${bed}.intersect_with_${feature_name}.piRNA.uniqueMappers.pp


Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features.R \
	${bed}.intersect_with_${feature_name}.unique_species \
	${bed}.intersect_with_${feature_name}.lendis \
	${bed}.intersect_with_${feature_name}.siRNA.lendis \
	${bed}.intersect_with_${feature_name}.piRNA.lendis \
	${ext_len} \
	${bed}.intersect_with_${feature_name}.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.uniqueMappers.pp \
	${bed}.intersect_with_${feature_name}.siRNA.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.siRNA.uniqueMappers.pp \
	${bed}.intersect_with_${feature_name}.piRNA.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.piRNA.uniqueMappers.pp \
1>&2 && \
Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features.R \
	${bed}.intersect_with_${feature_name}.all_reads \
	${bed}.intersect_with_${feature_name}.allMapper.lendis \
	${bed}.intersect_with_${feature_name}.allMapper.siRNA.lendis \
	${bed}.intersect_with_${feature_name}.allMapper.piRNA.lendis \
	${ext_len} \
	${bed}.intersect_with_${feature_name}.reads.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.pp \
	${bed}.intersect_with_${feature_name}.siRNA.reads.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.siRNA.pp \
	${bed}.intersect_with_${feature_name}.piRNA.reads.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.piRNA.pp \
1>&2 && \
Rscript $PIPELINE_DIRECTORY/bin/piPipes_draw_smallRNA_features_3p.R \
	${bed}.intersect_with_${feature_name}.3p \
	${bed}.intersect_with_${feature_name}.lendis \
	${bed}.intersect_with_${feature_name}.siRNA.lendis \
	${bed}.intersect_with_${feature_name}.piRNA.lendis \
	${ext_len} \
	${bed}.intersect_with_${feature_name}.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.species.3end_60.percentage \
	${bed}.intersect_with_${feature_name}.siRNA.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.siRNA.species.3end_60.percentage \
	${bed}.intersect_with_${feature_name}.piRNA.species.5end_60.percentage \
	${bed}.intersect_with_${feature_name}.piRNA.species.3end_60.percentage \
1>&2 && \
rm -rf ${bed}.intersect_with_${feature_name} ${bed}.intersect_with_${feature_name}.siRNA ${bed}.intersect_with_${feature_name}.piRNA
