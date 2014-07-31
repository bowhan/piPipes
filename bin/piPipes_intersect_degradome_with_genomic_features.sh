
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
GENOME_ALLMAP_BED12=$1 # genome wide alignment of all mappers, in bed12 format, with $5 being 1/NTM
SUMMARY=$2 # summary file
CPU=$3 # CPU to use
INTERSECT_OUTDIR=$4 # output directory to store files
SEED=${RANDOM}${RANDOM}${RANDOM}${RANDOM} # random name
. $COMMON_FOLDER/genomic_features # reading the information to intersect with, as well as some other annotation files
ALL_BED=`basename ${GENOME_ALLMAP_BED12%bed*}x_rpmk_MASK.bed12` # names for the file genernated here
# get rid of tRNA, rRNA, snoRNA...
bedtools_piPipes intersect -v -wa -split -a $GENOME_ALLMAP_BED12 -b $MASK | tee $INTERSECT_OUTDIR/${ALL_BED} | awk '{all_reads+=$5; if ($5==1) {unique_reads+=$4; ++unique_species}}END{printf "%d\t%d\t%d\n", unique_reads, all_reads, unique_species;}' > $INTERSECT_OUTDIR/.stats
print_header $SUMMARY
# doing intersecting and counting
para_file=$INTERSECT_OUTDIR/${SEED}.intersect.para
for t in ${TARGETS[@]}
do \
	echo "bash $DEBUG piPipes_degradome_intersect.sh $INTERSECT_OUTDIR/${ALL_BED} ${t} ${!t} $INTERSECT_OUTDIR/.stats $SRA_ALL_BED2" >> $para_file
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


