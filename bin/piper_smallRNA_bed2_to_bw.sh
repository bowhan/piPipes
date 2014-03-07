
# small RNA pipeline from pipipe: https://github.com/bowhan/pipipe.git
# pipipe: https://github.com/bowhan/pipipe.git
# An integrated pipeline for small RNA analysis 
# from small RNA Seq, RNASeq, CAGE/Degradome, ChIP-Seq and Genomic-Seq
# Wei Wang (wei.wang2@umassmed.edu)
# Bo W Han (bo.han@umassmed.edu, bowhan@me.com)
# the Zamore lab and the Weng lab
# University of Massachusetts Medical School

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
pipipe_bed2_to_bedGraph -i $OUTDIR/${INPUT_BED_NAME%bed2}sorted.bed2 -c $CHROM -p $CPU -o $OUTDIR/${INPUT_BED_NAME%bed2}sorted.
pipipe_bed2_to_bedGraph -i $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.bed2 -c $CHROM -p $CPU  -o $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.
# making bigWig
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Watson.bigWig" > $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph  > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.Crick.bigWig" >> $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Watson.bigWig" >> $para_file
echo -e "awk -v depth=$NormScale 'BEGIN{FS=OFS=\"\\\\t\"} {\$4*=depth; print}' $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph  > $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph.normalized && bedGraphToBigWig $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bedGraph.normalized $CHROM $OUTDIR/${INPUT_BED_NAME%bed2}sorted.uniq.Crick.bigWig" >> $para_file
ParaFly -c $para_file -CPU $CPU -failed_cmds ${para_file}.failedCommands 1>/dev/null 2>/dev/null && \
rm -rf $para_file ${para_file}".completed" # delete para file in case people rerun the pipeline