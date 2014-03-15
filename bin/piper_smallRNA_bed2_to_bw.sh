
# piper, a pipeline collection for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) <2014>  <Bo Han, Wei Wang, Phillip Zamore, Zhiping Weng>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


# this script takes 
INPUT_BED2=$1 # input bed2 format
INPUT_BED_NAME=`basename $INPUT_BED2`
CHROM=$2 # file with chrom size information
NormScale=$3 # normalization scale
CPU=$4 # CPU to use
OUTDIR=$5 && mkdir -p $OUTDIR # output dir
para_file=$OUTDIR/${RANDOM}${RANDOM}${RANDOM}.make_bigWig.para
# sorting is required for bedtools genomecov; using sort supporing MT
( sort --parallel=$CPU --temporary-directory=$OUTDIR -k1,1 $INPUT_BED2 | tee $OUTDIR/${INPUT_BED_NAME%bed2}sorted.bed2  | awk '$5==1' > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.bed2 ) || \
( sort                                               -k1,1 $INPUT_BED2 | tee $OUTDIR/${INPUT_BED_NAME%bed2}sorted.bed2  | awk '$5==1' > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.bed2 )
# making bedgraph
piper_bed2_to_bedGraph -i $OUTDIR/${INPUT_BED_NAME%bed2}sorted.bed2 -c $CHROM -p $CPU -o $OUTDIR/${INPUT_BED_NAME%bed2}sorted.
piper_bed2_to_bedGraph -i $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.bed2 -c $CHROM -p $CPU  -o $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.
# making bigWig
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bigWig" > $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph  > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bigWig" >> $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bigWig" >> $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph  > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bigWig" >> $para_file
ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>/dev/null 2>/dev/null && \
rm -rf $para_file ${para_file}".completed" # delete para file in case people rerun the pipeline